"""Protein–ligand contact screening using aligned receptor PDB and ligand SDF.

Standard amino-acid templates; no pH prediction, docking or energy inference.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
from .geometry import CartesianAtom, read_sdf
from .bond_orders import read_sdf_bonds
from .features import features_from_sdf
from .noncovalent import detect_non_covalent
from . import bond_distance, bond_angle, _add, _text
from xml.etree.ElementTree import Element, tostring

_AA=set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP".split())
_ACCEPT={"ASP":"OD1 OD2","GLU":"OE1 OE2","ASN":"OD1","GLN":"OE1","SER":"OG","THR":"OG1","TYR":"OH","MET":"SD","HID":"NE2","HIE":"ND1"}
_DONOR={"ARG":"NE NH1 NH2","LYS":"NZ","ASN":"ND2","GLN":"NE2","SER":"OG","THR":"OG1","TYR":"OH","TRP":"NE1","HID":"ND1","HIE":"NE2","HIP":"ND1 NE2"}
_HYD={"ALA":"CB","VAL":"CB CG1 CG2","LEU":"CB CG CD1 CD2","ILE":"CB CG1 CG2 CD1","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","TYR":"CB CG CD1 CD2 CE1 CE2","TRP":"CB CG CD1 CD2 CE2 CE3 CZ2 CZ3 CH2","PRO":"CB CG CD","MET":"CB CG CE"}
_RING={"PHE":[ "CG CD1 CE1 CZ CE2 CD2"],"TYR":["CG CD1 CE1 CZ CE2 CD2"],"TRP":["CG CD1 NE1 CE2 CD2","CD2 CE2 CZ2 CH2 CZ3 CE3"]}
for _h in ("HIS","HID","HIE","HIP"):
    _RING[_h]=["CG ND1 CE1 NE2 CD2"]
_METALS={"Zn","Mg","Ca","Fe","Mn","Cu","Co","Ni","Cd","Na","K"}

@dataclass(frozen=True)
class ReceptorAtom:
    atom: CartesianAtom
    residue: str
    name: str
    chain: str
    number: str
    insertion: str
    role: str

    @property
    def residue_id(self):
        return f"{self.chain or '_'}:{self.residue}{self.number}{self.insertion}"

@dataclass(frozen=True)
class ProteinContact:
    kind: str
    ligand_ids: tuple
    receptor_ids: tuple
    residue: str
    distance: float
    angle: float = None
    water_id: str = ""
    second_distance: float = None

def read_receptor(path):
    """First model, blank/A alternate locations; retain residue metadata."""
    result=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"): break
        if line[:6].strip() not in ("ATOM","HETATM") or line[16:17] not in (" ","A"): continue
        res=line[17:20].strip()
        element=line[76:78].strip().capitalize()
        if not element: raise ValueError("PDB element column required")
        role="protein" if res in _AA else "water" if res in ("HOH","WAT","DOD") else "metal" if element in _METALS else "other"
        atom=CartesianAtom("P:"+line[6:11].strip(),element,(line[30:38],line[38:46],line[46:54]))
        result.append(ReceptorAtom(atom,res,line[12:16].strip(),line[21:22].strip(),line[22:26].strip(),line[26:27].strip(),role))
    if not result or len({r.atom.id for r in result})!=len(result):
        raise ValueError("Empty receptor or duplicate atom serials")
    return result

def screen_protein_ligand(protein_pdb, ligand_sdf, *, record=0, pocket_cutoff=8.,
                          polar_cutoff=3.5, max_ha=2.5,min_angle=120., **thresholds):
    """Return geometric candidates; heavy-atom polar contacts are not H bonds.

    Standard side-chain protonation templates assume conventional charged
    Asp/Glu/Lys/Arg. HIS is left untyped for donor/acceptor roles unless HID/HIE/HIP.
    Water bridges without explicit H orientation are proximity candidates.
    """
    from math import isfinite
    if not all(isfinite(v) and v>0 for v in (pocket_cutoff,polar_cutoff,max_ha)) or not 0<=min_angle<=180:
        raise ValueError("Invalid screening thresholds")
    ligand=read_sdf(ligand_sdf,record=record)
    f=features_from_sdf(ligand_sdf,record=record)
    bonds=read_sdf_bonds(ligand_sdf,record=record)
    ligand=[CartesianAtom("L:"+a.id,a.element,a.xyz) for a in ligand]
    la={a.id:a for a in ligand}
    heavy=[a for a in ligand if a.element!="H"]
    if not heavy: raise ValueError("Ligand has no heavy atoms")
    raw=read_receptor(protein_pdb)
    nearby={r.residue_id for r in raw if r.role in ("protein","water","metal") and any(bond_distance(r.atom.xyz,a.xyz)<=pocket_cutoff for a in heavy)}
    rec=[r for r in raw if r.residue_id in nearby and r.role in ("protein","water","metal")]
    by={r.atom.id:r for r in rec}
    residues={}
    for r in rec:
        residues.setdefault(r.residue_id,{})[r.name]=r
    donor=set("L:"+i for i in f.donor_ids)
    accept=set("L:"+i for i in f.acceptor_ids)
    pd=set(); pa=set(); pos=set("L:"+i for i in f.positive_ids); neg=set("L:"+i for i in f.negative_ids)
    hyd=set("L:"+i for i in f.hydrophobic_ids)
    metal=set("L:"+i for i in f.metal_ids)
    coord=set("L:"+i for i in f.coordinating_ids)
    rings=[tuple("L:"+i for i in ring) for ring in f.aromatic_rings]
    for r in rec:
        i=r.atom.id; name=r.name; res=r.residue
        if r.role=="metal": metal.add(i)
        if r.role!="protein": continue
        if name in ("O","OXT") or name in _ACCEPT.get(res,"").split(): pa.add(i)
        if (name=="N" and res!="PRO") or name in _DONOR.get(res,"").split(): pd.add(i)
        if (res=="LYS" and name=="NZ") or (res=="ARG" and name in ("NH1","NH2")) or (res=="HIP" and name in ("ND1","NE2")): pos.add(i)
        if (res=="ASP" and name in ("OD1","OD2")) or (res=="GLU" and name in ("OE1","OE2")) or name=="OXT": neg.add(i)
        if name in _HYD.get(res,"").split(): hyd.add(i)
    coord.update(pa|{i for i in neg if i in by})
    for members in residues.values():
        res=next(iter(members.values())).residue
        for ring in _RING.get(res,[]):
            names=ring.split()
            if all(n in members for n in names): rings.append(tuple(members[n].atom.id for n in names))
    # Keep waters out of direct contact screening; group whole receptor as one partner.
    all_atoms=ligand+[r.atom for r in rec if r.role!="water"]
    groups={a.id:("ligand" if a.id.startswith("L:") else "receptor") for a in all_atoms}
    contacts=detect_non_covalent(all_atoms,positive_ids=pos,negative_ids=neg,hydrophobic_ids=hyd,
        aromatic_rings=rings,metal_ids=metal,coordinating_ids=coord,groups=groups,**thresholds)
    out=[]
    for c in contacts:
        l,p=(c.first,c.second) if c.first[0].startswith("L:") else (c.second,c.first)
        out.append(ProteinContact(c.kind,l,p,by[p[0]].residue_id,c.distance,c.angle))
    # Explicit donor H from ligand connectivity or nearest same-residue donor.
    xyz={a.id:a.xyz for a in ligand}; xyz.update({i:r.atom.xyz for i,r in by.items()})
    lh={}
    for b in bonds:
        a,c="L:"+b.start,"L:"+b.end
        if a in donor and la[c].element=="H": lh.setdefault(a,[]).append(c)
        if c in donor and la[a].element=="H": lh.setdefault(c,[]).append(a)
    ph={}
    for r in rec:
        if r.atom.element!="H" or r.role!="protein": continue
        candidates=[(bond_distance(r.atom.xyz,xyz[i]),i) for i in pd if by[i].residue_id==r.residue_id]
        if candidates:
            d,i=min(candidates)
            if .5<=d<=1.25: ph.setdefault(i,[]).append(r.atom.id)
    for ds,acs,hs,lig_donor in ((donor,pa,lh,True),(pd,accept,ph,False)):
        for d in sorted(ds):
            for a in sorted(acs):
                da=bond_distance(xyz[d],xyz[a])
                if not 0<da<=polar_cutoff: continue
                l,p=((d,),(a,)) if lig_donor else ((a,),(d,))
                found=False
                for h in hs.get(d,[]):
                    ha=bond_distance(xyz[h],xyz[a])
                    if 0<ha<=max_ha:
                        angle=bond_angle(xyz[d],xyz[h],xyz[a])
                        if angle>=min_angle:
                            out.append(ProteinContact("hydrogen",l,p,by[p[0]].residue_id,ha,angle))
                            found=True
                if not found and not hs.get(d):
                    out.append(ProteinContact("polar_proximity",l,p,by[p[0]].residue_id,da))
    for w in rec:
        if w.role!="water" or w.atom.element!="O": continue
        for l in sorted(donor|accept):
            dl=bond_distance(xyz[l],w.atom.xyz)
            if not 0<dl<=polar_cutoff: continue
            for p in sorted(pd|pa):
                dp=bond_distance(xyz[p],w.atom.xyz)
                if 0<dp<=polar_cutoff:
                    out.append(ProteinContact("water_bridge_candidate",(l,),(p,),by[p].residue_id,dl,
                                              water_id=w.atom.id,second_distance=dp))
    return sorted(out,key=lambda c:(c.kind,c.residue,c.ligand_ids,c.receptor_ids,c.distance))

def write_contact_csv(contacts,path):
    """Export atom identifiers, residue IDs and measured contact geometry."""
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(["kind","ligand_ids","receptor_ids","residue","distance_A","angle_deg","water_id","second_distance_A"])
        for c in contacts:
            writer.writerow([c.kind,";".join(c.ligand_ids),";".join(c.receptor_ids),c.residue,c.distance,c.angle,c.water_id,c.second_distance])
    return path

def protein_contact_svg(contacts, *, width=900,row_height=70,font_size=16):
    """Native SVG contact summary, not an atom-level molecular pose drawing."""
    contacts=list(contacts)
    if width<650 or row_height<55 or font_size<=0: raise ValueError("Invalid figure dimensions")
    height=80+row_height*max(1,len(contacts))
    root=Element("svg",{"xmlns":"http://www.w3.org/2000/svg","width":str(width),"height":str(height),"viewBox":f"0 0 {width} {height}"})
    _add(root,"rect",x=0,y=0,width=width,height=height,fill="white")
    _text(root,width/2,26,"Protein–ligand contact candidates","#17212E",font_size+4)
    if not contacts: _text(root,width/2,70,"No contacts detected","#17212E",font_size)
    for i,c in enumerate(contacts):
        y=60+i*row_height+row_height/2
        color={"hydrogen":"#16834E","salt_bridge":"#9333EA","hydrophobic":"#8A6914","metal_coordination":"#087E8B","water_bridge_candidate":"#2563EB"}.get(c.kind,"#526477")
        _text(root,width*.16,y,"Ligand "+",".join(c.ligand_ids),"#17212E",font_size)
        _text(root,width*.82,y,c.residue,"#17212E",font_size)
        _add(root,"line",x1=width*.32,y1=y,x2=width*.68,y2=y,stroke=color,stroke_width=2,stroke_dasharray="7 5")
        label=f"{c.distance:.2f} Å"
        if c.second_distance is not None: label+=f" / {c.second_distance:.2f} Å"
        _text(root,width*.5,y-15,c.kind.replace("_"," "),color,font_size)
        _text(root,width*.5,y+15,label,color,font_size)
    return tostring(root,encoding="unicode")

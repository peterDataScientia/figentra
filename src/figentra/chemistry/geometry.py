"""Read Cartesian geometries without inferring chemical connectivity."""
from dataclasses import dataclass
from pathlib import Path
from math import isfinite, hypot
from . import Atom, Bond, bond_distance, bond_angle, molecular_interaction

@dataclass(frozen=True)
class CartesianAtom:
    id: str
    element: str
    xyz: tuple

    def __post_init__(self):
        p=tuple(float(v) for v in self.xyz)
        if len(p)!=3 or not all(isfinite(v) for v in p):
            raise ValueError("Invalid Cartesian coordinates")
        object.__setattr__(self,"xyz",p)

def read_xyz(path):
    """Read a single XYZ frame; coordinates in angstroms."""
    lines=Path(path).read_text().splitlines()
    n=int(lines[0])
    if n<=0 or len(lines)<n+2 or any(s.strip() for s in lines[n+2:]):
        raise ValueError("Expected one complete XYZ frame")
    atoms=[]
    for i,line in enumerate(lines[2:n+2],1):
        fields=line.split()
        if len(fields)<4:
            raise ValueError("Invalid XYZ atom line")
        atoms.append(CartesianAtom(fields[0]+str(i),fields[0],tuple(fields[1:4])))
    return atoms

def read_pdb(path):
    """Read first PDB model, blank/A alternate locations, no bond inference."""
    atoms=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line[:6].strip() not in ("ATOM","HETATM") or line[16:17] not in (" ","A"):
            continue
        element=line[76:78].strip().capitalize()
        if not element:
            raise ValueError("PDB element column required")
        atoms.append(CartesianAtom(line[6:11].strip(),element,
                     (line[30:38],line[38:46],line[46:54])))
    if not atoms or len({a.id for a in atoms})!=len(atoms):
        raise ValueError("No atoms or duplicate serial numbers")
    return atoms

def read_sdf(path, *, record=0):
    """Read one SDF record via RDKit; preserve stored coordinates and hydrogens."""
    from rdkit import Chem
    supplier=Chem.SDMolSupplier(str(path),removeHs=False)
    if record<0 or record>=len(supplier):
        raise ValueError("SDF record outside range")
    mol=supplier[record]
    if mol is None or not mol.GetNumConformers():
        raise ValueError("Invalid SDF record or missing conformer")
    conf=mol.GetConformer()
    if not conf.Is3D():
        raise ValueError("SDF record is not marked as 3D")
    atoms=[CartesianAtom(str(a.GetIdx()+1),a.GetSymbol(),
            tuple(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms()]
    return atoms

def geometry_interaction(atoms, bonds, *, width=800,height=400,
                         atom_radius=20, projection="auto", annotate=True, **kwargs):
    """Project Cartesian atoms to a schematic and label true 3D distances.

    Bonds must be supplied explicitly. Choose xy/xz/yz or auto (maximise
    minimum projected separation). Crowded projections raise an error.
    """
    atoms,bonds=list(atoms),list(bonds)
    if not atoms or len({a.id for a in atoms})!=len(atoms):
        raise ValueError("Provide unique atoms")
    if not all(isfinite(v) and v>0 for v in (width,height,atom_radius)):
        raise ValueError("Invalid canvas dimensions")
    margin=atom_radius+45
    if min(width,height)<=2*margin:
        raise ValueError("Canvas too small")
    planes={"xy":(0,1),"xz":(0,2),"yz":(1,2)}
    def layout(plane):
        i,j=planes[plane]
        xs=[a.xyz[i] for a in atoms]; ys=[a.xyz[j] for a in atoms]
        dx=max(xs)-min(xs); dy=max(ys)-min(ys)
        scale=min((width-2*margin)/dx if dx else float("inf"),
                  (height-2*margin)/dy if dy else float("inf"))
        if not isfinite(scale): scale=1
        points=[(width/2+(x-(min(xs)+max(xs))/2)*scale,
                 height/2-(y-(min(ys)+max(ys))/2)*scale) for x,y in zip(xs,ys)]
        separation=min((hypot(x-u,y-v) for k,(x,y) in enumerate(points)
                        for u,v in points[k+1:]),default=float("inf"))
        return separation,points
    if projection not in ("auto",*planes):
        raise ValueError("Invalid projection")
    separation,points=max((layout(p) for p in planes),key=lambda r:r[0]) if projection=="auto" else layout(projection)
    if separation<=2*atom_radius:
        raise ValueError("Crowded projection: select fewer atoms or a larger canvas")
    by_id={a.id:a for a in atoms}
    drawn=[Atom(a.id,a.element,x,y,a.element) for a,(x,y) in zip(atoms,points)]
    annotated=[]
    for b in bonds:
        if b.start not in by_id or b.end not in by_id:
            raise ValueError("Unknown bond endpoint")
        label=b.label or (f"{bond_distance(by_id[b.start].xyz,by_id[b.end].xyz):.3f} Å" if annotate else "")
        annotated.append(Bond(b.start,b.end,b.kind,label))
    return molecular_interaction(drawn,annotated,width=width,height=height,
                                 atom_radius=atom_radius,**kwargs)

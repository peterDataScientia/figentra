"""Native SVG central-ligand contact maps using RDKit 2D depiction.

Layout coordinates are schematic. Reported distances are retained from the
3D contact calculation, never measured from the drawing.
"""
from math import isfinite
from xml.etree import ElementTree as ET
from . import _add, _text

_COLORS={"hydrogen":"#16834E","polar_proximity":"#526477",
         "salt_bridge":"#9333EA","hydrophobic":"#8A6914",
         "metal_coordination":"#087E8B","pi_parallel":"#B45309",
         "pi_t_shaped":"#B45309","cation_pi":"#9F357A",
         "water_bridge_candidate":"#2563EB"}

def protein_ligand_map(ligand_sdf, contacts, *, record=0, width=1200,
                       height=800,font_size=16,max_contacts=24,
                       background="white", auto_height=True, show_legend=True):
    """Draw a central ligand and contact cards in two non-overlapping columns.

    Each contact is retained separately; if max_contacts is exceeded the
    caller must explicitly select a subset or enlarge the limit. Ring contacts
    anchor to the mean of their member atom depiction coordinates. Water nodes
    are schematic and their two distances are kept separate.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    contacts=list(contacts)
    if not all(isfinite(v) and v>0 for v in (width,height,font_size)) or width<900 or height<500 or font_size>24:
        raise ValueError("Use width >=900, height >=500 and font size <=24")
    if len(contacts)>max_contacts:
        raise ValueError("Too many contacts: select a subset or increase max_contacts")
    if background != "white":
        raise ValueError("Contact map currently supports white background for contrast")
    cardh=max(66,font_size*4.2)
    gap=14
    legend_rows=(len({c.kind for c in contacts})+2)//3 if show_legend else 0
    footer=65+legend_rows*30
    required=160+footer+((len(contacts)+1)//2)*(cardh+gap)
    if auto_height:
        height=max(height,required)
    supplier=Chem.SDMolSupplier(str(ligand_sdf),removeHs=False)
    if record<0 or record>=len(supplier) or supplier[record] is None:
        raise ValueError("Invalid SDF record")
    original=supplier[record]
    for a in original.GetAtoms():
        a.SetIntProp("_figentra_index",a.GetIdx()+1)
    mol=Chem.RemoveHs(original)
    rdDepictor.Compute2DCoords(mol)
    iw=int(width*.40); ih=int(height*.58)
    ox=(width-iw)/2; oy=(height-ih)/2
    drawer=rdMolDraw2D.MolDraw2DSVG(iw,ih)
    opts=drawer.drawOptions()
    opts.clearBackground=False
    opts.bondLineWidth=2.5
    opts.padding=.12
    opts.addStereoAnnotation=True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    points={}
    for a in mol.GetAtoms():
        p=drawer.GetDrawCoords(a.GetIdx())
        points["L:"+str(a.GetIntProp("_figentra_index"))]=(p.x+ox,p.y+oy)
    root=ET.Element("svg",{"xmlns":"http://www.w3.org/2000/svg",
        "width":str(width),"height":str(height),"viewBox":f"0 0 {width} {height}"})
    _add(root,"rect",x=0,y=0,width=width,height=height,fill=background)
    _text(root,width/2,30,"Protein–ligand interaction map","#17212E",font_size+6)
    _text(root,width/2,57,"2D schematic • distances from 3D pose","#526477",font_size-2,"normal")
    sides=[[],[]]
    for c in contacts:
        if not c.ligand_ids or any(i not in points for i in c.ligand_ids):
            raise ValueError("Contact refers to missing or hidden ligand atom: "+str(c.ligand_ids))
        if c.kind not in _COLORS:
            raise ValueError("Unsupported contact kind: "+c.kind)
        if not isfinite(c.distance) or c.distance<=0:
            raise ValueError("Invalid contact distance")
        if c.kind=="water_bridge_candidate" and (not c.water_id or c.second_distance is None or not isfinite(c.second_distance) or c.second_distance<=0):
            raise ValueError("Water bridge requires both distances and water ID")
        x=sum(points[i][0] for i in c.ligand_ids)/len(c.ligand_ids)
        y=sum(points[i][1] for i in c.ligand_ids)/len(c.ligand_ids)
        sides[0 if x<width/2 else 1].append((y,x,c))
    # Balance column counts by moving contacts nearest the central axis.
    while abs(len(sides[0])-len(sides[1]))>1:
        source=0 if len(sides[0])>len(sides[1]) else 1
        k=min(range(len(sides[source])),key=lambda j:abs(sides[source][j][1]-width/2))
        sides[1-source].append(sides[source].pop(k))
    cardw=width*.23
    available=height-160-footer
    if any(len(s)*cardh+max(0,len(s)-1)*gap>available for s in sides):
        raise ValueError("Crowded contact map: increase height or select contacts")
    # Place connectors first, then the chemical drawing and cards.
    cards=[]
    water_nodes=[]
    for side,items in enumerate(sides):
        items.sort(key=lambda p:(p[0],p[1],p[2].residue,p[2].kind))
        total=len(items)*cardh+max(0,len(items)-1)*gap
        for j,(ay,ax,c) in enumerate(items):
            cy=(height-footer-total)/2+cardh/2+j*(cardh+gap)+30
            cx=cardw/2+20 if side==0 else width-cardw/2-20
            edge=cx+cardw/2 if side==0 else cx-cardw/2
            color=_COLORS[c.kind]
            rail=ox-18 if side==0 else ox+iw+18
            # Fan out toward a routing rail outside the ligand drawing.
            _add(root,"path",d=f"M{ax} {ay} L{rail} {cy}",fill="none",stroke=color,stroke_width=2,stroke_dasharray="6 5")
            if c.kind=="water_bridge_candidate":
                wx=(rail+edge)/2; wy=cy
                _add(root,"line",x1=rail,y1=cy,x2=wx,y2=wy,stroke=color,stroke_width=2,stroke_dasharray="6 5")
                _add(root,"line",x1=wx,y1=wy,x2=edge,y2=cy,stroke=color,stroke_width=2,stroke_dasharray="6 5")
                water_nodes.append((wx,wy,color))
            else:
                _add(root,"line",x1=rail,y1=cy,x2=edge,y2=cy,stroke=color,stroke_width=2,stroke_dasharray="6 5")
            cards.append((cx,cy,c,color))
    ligand_svg=ET.fromstring(drawer.GetDrawingText())
    ligand_svg.set("x",str(ox)); ligand_svg.set("y",str(oy))
    ligand_svg.set("width",str(iw)); ligand_svg.set("height",str(ih))
    root.append(ligand_svg)
    for wx,wy,color in water_nodes:
        _add(root,"circle",cx=wx,cy=wy,r=17,fill="white",stroke=color,stroke_width=2)
        _text(root,wx,wy,"H₂O",color,12)
    for cx,cy,c,color in cards:
        _add(root,"rect",x=cx-cardw/2,y=cy-cardh/2,width=cardw,height=cardh,rx=8,fill="white",stroke=color,stroke_width=1.8)
        _text(root,cx,cy-font_size*1.25,c.residue,"#17212E",font_size)
        _text(root,cx,cy,c.kind.replace("_"," "),color,font_size*.82)
        distance=f"{c.distance:.2f} Å"
        if c.kind=="hydrogen": distance+=" (H···A)"
        if c.second_distance is not None: distance+=f" / {c.second_distance:.2f} Å"
        _text(root,cx,cy+font_size*1.25,distance,color,font_size*.85)
    if show_legend:
        for j,kind in enumerate(sorted({c.kind for c in contacts})):
            x=width*(.17+(j%3)*.33); y=height-footer+25+(j//3)*30
            _add(root,"line",x1=x-width*.13,y1=y,x2=x-width*.10,y2=y,stroke=_COLORS[kind],stroke_width=2,stroke_dasharray="6 5")
            _text(root,x+width*.02,y,kind.replace("_"," "),_COLORS[kind],font_size*.85)
    if not contacts:
        _text(root,width/2,height-70,"No contacts detected","#526477",font_size)
    _text(root,width/2,height-32,"Water bridges and proximity contacts are geometric candidates.","#526477",font_size-2,"normal")
    return ET.tostring(root,encoding="unicode")

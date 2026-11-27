"""Native SVG molecular schematics. Coordinates are illustrative unless supplied as geometry."""
from dataclasses import dataclass
from math import dist, acos, degrees, isfinite, hypot
from pathlib import Path
from xml.etree.ElementTree import Element, SubElement, tostring

__all__ = ["Atom", "Bond", "molecular_interaction", "hydrogen_bond_network",
           "bond_distance", "bond_angle", "save_svg"]

@dataclass(frozen=True)
class Atom:
    id: str
    element: str
    x: float
    y: float
    label: str = ""

@dataclass(frozen=True)
class Bond:
    start: str
    end: str
    kind: str = "single"
    label: str = ""

_COLORS = {"O": "#D9232B", "N": "#225ACB", "C": "#344054",
           "H": "#E5E7EB", "S": "#D5A21C", "P": "#D97706"}
_KINDS = {"single", "double", "triple", "aromatic", "hydrogen", "forming", "breaking", "salt_bridge", "hydrophobic", "metal_coordination"}

def _point(p):
    p = tuple(float(x) for x in p)
    if len(p) != 3 or not all(isfinite(x) for x in p):
        raise ValueError("Expected three finite Cartesian coordinates")
    return p

def bond_distance(a, b):
    """Distance in the same units as the input Cartesian coordinates."""
    return dist(_point(a), _point(b))

def bond_angle(a, b, c):
    """Angle ABC in degrees, using three-dimensional coordinates."""
    a,b,c = map(_point, (a,b,c))
    u = [x-y for x,y in zip(a,b)]
    v = [x-y for x,y in zip(c,b)]
    norm = hypot(*u)*hypot(*v)
    if not norm:
        raise ValueError("Angle undefined for coincident points")
    return degrees(acos(max(-1., min(1., sum(x*y for x,y in zip(u,v))/norm))))

def _add(parent, tag, **attrs):
    return SubElement(parent, tag, {k.replace("_","-"):str(v) for k,v in attrs.items()})

def _text(root, x, y, value, color, size, weight="bold"):
    el = _add(root, "text", x=x, y=y, fill=color, font_size=size,
              font_family="DejaVu Sans, sans-serif", font_weight=weight,
              text_anchor="middle", dominant_baseline="central")
    el.text = str(value)

def molecular_interaction(atoms, bonds, *, width=800, height=360,
                          background="white", text_color=None, font_size=18,
                          atom_radius=25, line_width=4, caption=""):
    """Return standalone SVG. Atom x/y are layout pixels, not physical coordinates.

    Bond labels are caller-supplied; no connectivity or chemistry is inferred.
    """
    atoms, bonds = list(atoms), list(bonds)
    if not atoms or len({a.id for a in atoms}) != len(atoms):
        raise ValueError("Provide atoms with unique IDs")
    if not all(isfinite(v) and v>0 for v in (width,height,font_size,atom_radius,line_width)):
        raise ValueError("Canvas and style dimensions must be positive and finite")
    by_id = {a.id:a for a in atoms}
    radius = {a.id:atom_radius*(.70 if a.element=="H" else 1) for a in atoms}
    for a in atoms:
        r=radius[a.id]
        if not (r <= a.x <= width-r and r <= a.y <= height-r):
            raise ValueError("Atom outside canvas: "+a.id)
    for b in bonds:
        if b.kind not in _KINDS or b.start not in by_id or b.end not in by_id or b.start==b.end:
            raise ValueError("Invalid bond: "+str(b))
    dark = background.lower() in ("black", "#202020", "#222222", "#111111")
    fg = text_color or ("#FFFFFF" if dark else "#17212E")
    root = Element("svg", {"xmlns":"http://www.w3.org/2000/svg", "width":str(width),
                           "height":str(height),"viewBox":f"0 0 {width} {height}"})
    _add(root,"rect",x=0,y=0,width=width,height=height,fill=background)
    for b in bonds:
        a,c=by_id[b.start],by_id[b.end]
        dx,dy=c.x-a.x,c.y-a.y
        d=hypot(dx,dy)
        if d <= radius[a.id]+radius[c.id]:
            raise ValueError("Bonded atoms overlap: "+str(b))
        ux,uy=dx/d,dy/d
        x1,y1=a.x+ux*radius[a.id],a.y+uy*radius[a.id]
        x2,y2=c.x-ux*radius[c.id],c.y-uy*radius[c.id]
        color={"hydrogen":"#16834E","forming":"#2563EB","breaking":"#B45309","salt_bridge":"#9333EA","hydrophobic":"#8A6914","metal_coordination":"#087E8B"}.get(b.kind,"#9299A3")
        offsets={"double":(-3,3),"triple":(-6,0,6),"aromatic":(-3,3)}.get(b.kind,(0,))
        for o in offsets:
            attrs=dict(x1=x1-uy*o,y1=y1+ux*o,x2=x2-uy*o,y2=y2+ux*o,
                       stroke=color,stroke_width=line_width,stroke_linecap="round")
            if b.kind in ("hydrogen","forming","breaking","salt_bridge","hydrophobic","metal_coordination") or (b.kind=="aromatic" and o>0):
                attrs["stroke_dasharray"]="8 7"
            _add(root,"line",**attrs)
        if b.label:
            nx,ny=-uy,ux
            if ny>0: nx,ny=-nx,-ny
            _text(root,(x1+x2)/2+nx*20,(y1+y2)/2+ny*20,b.label,
                  color if b.kind in ("hydrogen","forming","breaking","salt_bridge","hydrophobic","metal_coordination") else fg,font_size)
    for a in atoms:
        _add(root,"circle",cx=a.x,cy=a.y,r=radius[a.id],
             fill=_COLORS.get(a.element,"#526477"),stroke="#87909B",stroke_width=1)
        _text(root,a.x,a.y,a.label or a.id,"#222222" if a.element=="H" else "#FFFFFF",font_size)
    if caption:
        _text(root,width/2,height-24,caption,fg,font_size,"normal")
    return tostring(root,encoding="unicode")

def hydrogen_bond_network(atoms, bonds, **kwargs):
    """Draw explicitly specified hydrogen bonds and covalent connections."""
    bonds=list(bonds)
    if not any(b.kind=="hydrogen" for b in bonds):
        raise ValueError("Specify at least one hydrogen bond")
    return molecular_interaction(atoms,bonds,**kwargs)

def save_svg(svg, path):
    """Save native editable SVG source."""
    path=Path(path)
    if path.suffix.lower() != ".svg":
        raise ValueError("Expected .svg extension")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(svg,encoding="utf-8")
    return path

from .reactions import (ElectronArrow, reaction_mechanism, molecular_structure,
                        reaction_scheme, protonation_states, export_chemistry)
__all__ += ["ElectronArrow", "reaction_mechanism", "molecular_structure",
            "reaction_scheme", "protonation_states", "export_chemistry"]

from .geometry import CartesianAtom, read_xyz, read_pdb, read_sdf, geometry_interaction
__all__ += ["CartesianAtom", "read_xyz", "read_pdb", "read_sdf", "geometry_interaction"]

from .detection import detect_bonds, detect_covalent_bonds, detect_hydrogen_bonds
__all__ += ["detect_bonds", "detect_covalent_bonds", "detect_hydrogen_bonds"]

from .bond_orders import assign_bond_orders, read_sdf_bonds
__all__ += ["assign_bond_orders", "read_sdf_bonds"]

from .noncovalent import Interaction, detect_non_covalent, interaction_bonds
__all__ += ["Interaction", "detect_non_covalent", "interaction_bonds"]

from .features import ChemicalFeatures, perceive_features, features_from_smiles, features_from_sdf, screen_sdf_contacts
__all__ += ["ChemicalFeatures", "perceive_features", "features_from_smiles", "features_from_sdf", "screen_sdf_contacts"]

from .protein import ReceptorAtom, ProteinContact, read_receptor, screen_protein_ligand, write_contact_csv, protein_contact_svg
__all__ += ["ReceptorAtom", "ProteinContact", "read_receptor", "screen_protein_ligand", "write_contact_csv", "protein_contact_svg"]

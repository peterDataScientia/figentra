"""SVG reaction schemes and optional RDKit/CairoSVG adapters."""
from pathlib import Path
from xml.etree import ElementTree as ET
from dataclasses import dataclass
from math import isfinite
from . import _add, _text, save_svg

_NS="http://www.w3.org/2000/svg"
ET.register_namespace("", _NS)

@dataclass(frozen=True)
class ElectronArrow:
    start: tuple
    control: tuple
    end: tuple
    label: str = ""
    color: str = "#2563EB"

def reaction_mechanism(svg, arrows):
    """Overlay caller-defined two-electron curved arrows on a native SVG.
    Coordinates refer to the SVG viewBox. No mechanism is inferred.
    """
    root=ET.fromstring(svg)
    if root.tag not in ("svg", "{"+_NS+"}svg"):
        raise ValueError("Expected SVG document")
    def add(parent, tag, **attrs):
        return ET.SubElement(parent, "{"+_NS+"}"+tag,
                              {k.replace("_","-"):str(v) for k,v in attrs.items()})
    defs=add(root,"defs")
    for i,a in enumerate(arrows):
        if any(len(p)!=2 or not all(isfinite(v) for v in p) for p in (a.start,a.control,a.end)):
            raise ValueError("Arrow requires three finite 2D points")
        key=f"figentra-electron-arrow-{i}"
        marker=add(defs,"marker",id=key,marker_width=7,marker_height=7,
                   refX=6,refY=3.5,orient="auto",viewBox="0 0 7 7")
        add(marker,"path",d="M0 0 L7 3.5 L0 7 Z",fill=a.color)
        x,y=a.start; cx,cy=a.control; ex,ey=a.end
        add(root,"path",d=f"M{x} {y} Q{cx} {cy} {ex} {ey}",
            fill="none",stroke=a.color,stroke_width=2.5,marker_end=f"url(#{key})")
        if a.label:
            t=add(root,"text",x=cx,y=cy-10,fill=a.color,font_size=16,
                  font_weight="bold",text_anchor="middle")
            t.text=a.label
    return ET.tostring(root,encoding="unicode")

def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise ImportError("Install RDKit to draw chemical structures") from exc
    return Chem,rdMolDraw2D

def molecular_structure(smiles, *, width=400, height=300):
    """Generate a chemically parsed 2D structure from SMILES using RDKit."""
    Chem,draw=_rdkit()
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    drawer=draw.MolDraw2DSVG(width,height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def reaction_scheme(reaction_smiles, *, width=1000,height=300):
    """Draw explicitly supplied reactants>agents>products; no prediction."""
    Chem,draw=_rdkit()
    from rdkit.Chem import rdChemReactions
    rxn=rdChemReactions.ReactionFromSmarts(reaction_smiles,useSmiles=True)
    if rxn is None or not rxn.GetNumReactantTemplates() or not rxn.GetNumProductTemplates():
        raise ValueError("Expected reaction SMILES with reactants and products")
    drawer=draw.MolDraw2DSVG(width,height)
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def protonation_states(states, *, cell_width=300,height=280):
    """Compare named, explicitly provided SMILES states; no pKa inference."""
    Chem,draw=_rdkit()
    states=list(states)
    if not states:
        raise ValueError("Provide (label, SMILES) pairs")
    mols=[Chem.MolFromSmiles(s) for _,s in states]
    if any(m is None for m in mols):
        raise ValueError("Invalid state SMILES")
    drawer=draw.MolDraw2DSVG(cell_width*len(states),height,cell_width,height)
    drawer.DrawMolecules(mols,legends=[str(label) for label,_ in states])
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def export_chemistry(svg,path,*,dpi=600):
    """Export SVG directly or render PDF/PNG/TIFF with optional dependencies.

    SVG dimensions are interpreted as CSS pixels (96 pixels per inch).
    """
    path=Path(path)
    ext=path.suffix.lower()
    if ext==".svg":
        return save_svg(svg,path)
    if ext not in (".pdf",".png",".tif",".tiff"):
        raise ValueError("Expected SVG, PDF, PNG or TIFF")
    if not isfinite(dpi) or dpi<=0:
        raise ValueError("DPI must be positive")
    try:
        import cairosvg
    except ImportError as exc:
        raise ImportError("Install CairoSVG and system Cairo for conversion") from exc
    path.parent.mkdir(parents=True,exist_ok=True)
    data=svg.encode("utf-8")
    if ext==".pdf":
        cairosvg.svg2pdf(bytestring=data,write_to=str(path))
    else:
        from io import BytesIO
        from PIL import Image
        png=cairosvg.svg2png(bytestring=data,scale=dpi/96)
        with Image.open(BytesIO(png)) as im:
            im.save(path,dpi=(dpi,dpi))
    return path

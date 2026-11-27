import math
import xml.etree.ElementTree as ET
import pytest
from figentra.chemistry import (
    Atom, Bond, ElectronArrow, molecular_interaction,
    hydrogen_bond_network, bond_distance, bond_angle,
    reaction_mechanism, molecular_structure, reaction_scheme,
    protonation_states, export_chemistry,
)

NS={"s":"http://www.w3.org/2000/svg"}

def example():
    return hydrogen_bond_network(
        [Atom("O1","O",80,100), Atom("H1","H",180,100),
         Atom("O2","O",320,100)],
        [Bond("O1","H1","single","0.988 Å"),
         Bond("H1","O2","hydrogen","1.775 Å")],
        width=400,height=220)

def test_geometry():
    assert bond_distance((0,0,0),(1,2,2)) == pytest.approx(3)
    assert bond_angle((1,0,0),(0,0,0),(0,1,0)) == pytest.approx(90)
    assert bond_angle((-1,0,0),(0,0,0),(1,0,0)) == pytest.approx(180)
    with pytest.raises(ValueError):
        bond_angle((0,0,0),(0,0,0),(1,0,0))

def test_native_svg():
    root=ET.fromstring(example())
    assert len(root.findall("s:circle",NS))==3
    assert len(root.findall("s:line",NS))==2
    assert any(x.get("stroke-dasharray") for x in root.findall("s:line",NS))
    assert "1.775 Å" in "".join(root.itertext())

def test_arrow():
    root=ET.fromstring(reaction_mechanism(example(),[
        ElectronArrow((80,60),(150,10),(240,60))]))
    assert root.find("s:defs/s:marker",NS) is not None
    assert any(x.get("marker-end") for x in root.findall("s:path",NS))

def test_svg_export(tmp_path):
    p=export_chemistry(example(),tmp_path/"figure.svg")
    assert ET.parse(p).getroot().tag.endswith("svg")

def test_rdkit():
    pytest.importorskip("rdkit")
    for svg in (molecular_structure("CCO"),
                reaction_scheme("CCO>>CC=O"),
                protonation_states([("neutral","CN"),("cation","C[NH3+]")])):
        assert ET.fromstring(svg).tag.endswith("svg")
    with pytest.raises(ValueError):
        molecular_structure("invalid smiles")

@pytest.mark.parametrize("ext",["png","pdf","tiff"])
def test_conversion(tmp_path,ext):
    pytest.importorskip("cairosvg")
    p=export_chemistry(example(),tmp_path/("figure."+ext),dpi=96)
    assert p.stat().st_size>0

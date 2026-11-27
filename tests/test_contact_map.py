import xml.etree.ElementTree as ET
import pytest
from figentra.chemistry import ProteinContact, protein_ligand_map

@pytest.fixture
def ligand(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    mol=Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    rdDepictor.Compute2DCoords(mol)
    path=tmp_path/"ligand.sdf"
    w=Chem.SDWriter(str(path)); w.write(mol); w.close()
    return path

def test_map(ligand):
    contacts=[
        ProteinContact("hydrogen",("L:7",),("P:1",),"A:ASN140",1.95,165.2),
        ProteinContact("pi_parallel",tuple("L:"+str(i) for i in range(1,7)),("P:2",),"A:TYR97",4.23),
        ProteinContact("water_bridge_candidate",("L:7",),("P:3",),"A:GLN85",2.81,water_id="P:4",second_distance=2.94),
    ]
    root=ET.fromstring(protein_ligand_map(ligand,contacts))
    text=" ".join(root.itertext())
    for s in ("A:ASN140","A:TYR97","A:GLN85","1.95","2.81","2.94","H₂O"):
        assert s in text

def test_invalid_anchor(ligand):
    c=ProteinContact("hydrogen",("L:999",),("P:1",),"A:ASN140",2.0)
    with pytest.raises(ValueError,match="missing or hidden"):
        protein_ligand_map(ligand,[c])

def test_empty(ligand):
    assert "No contacts detected" in protein_ligand_map(ligand,[])

def test_crowding(ligand):
    c=ProteinContact("hydrogen",("L:7",),("P:1",),"A:ASN140",2.0)
    with pytest.raises(ValueError,match="Crowded"):
        protein_ligand_map(ligand,[c]*12,height=500,auto_height=False)


def test_auto_height(ligand):
    c=ProteinContact("hydrogen",("L:7",),("P:1",),"A:ASN140",2.0)
    root=ET.fromstring(protein_ligand_map(ligand,[c]*12,height=500))
    assert float(root.get("height"))>500

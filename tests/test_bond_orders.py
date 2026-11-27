import pytest
import xml.etree.ElementTree as ET
from figentra.chemistry import CartesianAtom, assign_bond_orders, read_sdf_bonds, Atom, Bond, molecular_interaction

@pytest.mark.parametrize("smiles,expected",[("CC","single"),("C=C","double"),("C#C","triple"),("c1ccccc1","aromatic"),("[NH4+]","single")])
def test_inference(smiles,expected):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol,randomSeed=42)==0
    conf=mol.GetConformer()
    atoms=[CartesianAtom(str(a.GetIdx()),a.GetSymbol(),tuple(conf.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms()]
    bonds=assign_bond_orders(atoms,charge=Chem.GetFormalCharge(mol))
    assert len(bonds)==mol.GetNumBonds()
    assert expected in {b.kind for b in bonds}

def test_aromatic_svg():
    svg=molecular_interaction([Atom("C1","C",80,100),Atom("C2","C",220,100)],[Bond("C1","C2","aromatic")],width=300,height=200)
    root=ET.fromstring(svg)
    lines=root.findall("{http://www.w3.org/2000/svg}line")
    assert len(lines)==2
    assert sum(bool(l.get("stroke-dasharray")) for l in lines)==1

def test_sdf_orders(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    mol=Chem.MolFromSmiles("c1ccccc1")
    path=tmp_path/"benzene.sdf"
    writer=Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    bonds=read_sdf_bonds(path)
    assert len(bonds)==6
    assert all(b.kind=="aromatic" for b in bonds)

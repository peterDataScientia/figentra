import xml.etree.ElementTree as ET
import pytest
from figentra.chemistry import (
    CartesianAtom, Bond, read_xyz, read_pdb, read_sdf,
    geometry_interaction, bond_distance,
)

def test_xyz_and_projection(tmp_path):
    path=tmp_path/"water.xyz"
    path.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    atoms=read_xyz(path)
    assert len(atoms)==3
    assert bond_distance(atoms[0].xyz,atoms[1].xyz)==pytest.approx(.96)
    bonds=[Bond("O1","H2"),Bond("O1","H3")]
    for width in (600,1000):
        svg=geometry_interaction(atoms,bonds,width=width,height=450)
        assert "0.960 Å" in "".join(ET.fromstring(svg).itertext())

def test_pdb(tmp_path):
    path=tmp_path/"water.pdb"
    path.write_text("HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 20.00           O  \nEND\n")
    a=read_pdb(path)[0]
    assert a.element=="O"
    assert a.xyz==(1.,2.,3.)

def test_sdf(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles("O"))
    assert AllChem.EmbedMolecule(mol,randomSeed=42)==0
    path=tmp_path/"water.sdf"
    writer=Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    assert len(read_sdf(path))==3

def test_crowded_projection():
    atoms=[CartesianAtom("a","O",(0,0,0)),CartesianAtom("b","H",(0,0,0))]
    with pytest.raises(ValueError,match="Crowded"):
        geometry_interaction(atoms,[Bond("a","b")])

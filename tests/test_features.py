import pytest
pytest.importorskip("rdkit")
from rdkit import Chem
from figentra.chemistry import features_from_smiles, perceive_features, features_from_sdf

def test_amide():
    f=features_from_smiles("CC(=O)N")
    assert "3" in f.acceptor_ids
    assert "4" in f.donor_ids
    assert "4" not in f.acceptor_ids

def test_formal_charges():
    f=features_from_smiles("C[NH3+].CC(=O)[O-]")
    assert f.positive_ids==("2",)
    assert f.negative_ids==("6",)

def test_aromatic_and_halogen():
    f=features_from_smiles("Clc1ccccc1")
    assert len(f.aromatic_rings)==1
    assert len(f.aromatic_rings[0])==6
    assert f.halogen_donors==(("2","1"),)

def test_custom_ids():
    f=features_from_smiles("CC(=O)N",atom_ids=["a","b","c","d"])
    assert "c" in f.acceptor_ids
    with pytest.raises(ValueError):
        features_from_smiles("CC",atom_ids=["a","a"])

def test_sdf_atom_mapping(tmp_path):
    mol=Chem.AddHs(Chem.MolFromSmiles("CC(=O)N"))
    path=tmp_path/"compound.sdf"
    writer=Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    assert features_from_sdf(path)==perceive_features(mol)

def test_invalid_smiles():
    with pytest.raises(ValueError):
        features_from_smiles("invalid smiles")

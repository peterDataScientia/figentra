import pytest
from figentra.chemistry import CartesianAtom, Bond, detect_bonds, detect_covalent_bonds, detect_hydrogen_bonds

def test_water_dimer():
    atoms=[CartesianAtom("O1","O",(0,0,0)),CartesianAtom("H2","H",(.96,0,0)),
           CartesianAtom("H3","H",(-.24,.93,0)),CartesianAtom("O4","O",(2.75,0,0)),
           CartesianAtom("H5","H",(3.71,0,0)),CartesianAtom("H6","H",(2.51,.93,0))]
    cov=detect_covalent_bonds(atoms)
    assert len(cov)==4
    hb=detect_hydrogen_bonds(atoms,cov,acceptor_ids=["O1","O4"])
    assert [(b.start,b.end) for b in hb]==[("H2","O4")]
    assert hb[0].label=="1.790 Å"
    assert len(detect_bonds(atoms,acceptor_ids=["O1","O4"]))==5

def test_unsupported_element():
    with pytest.raises(ValueError,match="Unsupported"):
        detect_covalent_bonds([CartesianAtom("Fe1","Fe",(0,0,0))])

def test_angle_filter():
    atoms=[CartesianAtom("O","O",(0,0,0)),CartesianAtom("H","H",(1,0,0)),
           CartesianAtom("A","O",(1,2,0))]
    assert detect_hydrogen_bonds(atoms,[Bond("O","H")],acceptor_ids=["A"])==[]

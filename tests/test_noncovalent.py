import math
import pytest
from figentra.chemistry import CartesianAtom, detect_non_covalent, interaction_bonds

def atom(i,e,x,y,z):
    return CartesianAtom(i,e,(x,y,z))

def test_atom_contacts():
    atoms=[atom("p","N",0,0,0),atom("n","O",3,0,0),
           atom("c1","C",0,10,0),atom("c2","C",3.5,10,0),
           atom("zn","Zn",0,20,0),atom("o","O",2.1,20,0)]
    contacts=detect_non_covalent(atoms,positive_ids=["p"],negative_ids=["n"],
        hydrophobic_ids=["c1","c2"],metal_ids=["zn"],coordinating_ids=["o"])
    assert {c.kind for c in contacts}=={"salt_bridge","hydrophobic","metal_coordination"}
    assert len(interaction_bonds(contacts))==3

def test_exclusions():
    atoms=[atom("p","N",0,0,0),atom("n","O",3,0,0)]
    assert not detect_non_covalent(atoms,positive_ids=["p"],negative_ids=["n"],
        excluded_pairs=[("p","n")])
    assert not detect_non_covalent(atoms,positive_ids=["p"],negative_ids=["n"],
        groups={"p":"ligand","n":"ligand"})

def test_parallel_rings_and_cation():
    atoms=[]; rings=[]
    for k,z in enumerate((0,3.5)):
        ring=[]
        for j in range(6):
            i=f"r{k}_{j}"; t=j*math.pi/3
            atoms.append(atom(i,"C",1.4*math.cos(t),1.4*math.sin(t),z)); ring.append(i)
        rings.append(ring)
    atoms.append(atom("p","N",0,0,-3.5))
    contacts=detect_non_covalent(atoms,aromatic_rings=rings,positive_ids=["p"])
    assert {c.kind for c in contacts}=={"pi_parallel","cation_pi"}
    assert all(c.distance==pytest.approx(3.5) for c in contacts)
    with pytest.raises(ValueError,match="centroid"):
        interaction_bonds(contacts)

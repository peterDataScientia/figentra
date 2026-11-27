import csv
import xml.etree.ElementTree as ET
import pytest
from figentra.chemistry import read_receptor, ProteinContact, write_contact_csv, protein_contact_svg, screen_protein_ligand

def pdb_line(serial,name,res,chain,num,x,y,z,element,record="ATOM"):
    return f"{record:<6}{serial:5d} {name:>4s} {res:>3s} {chain}{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}\n"

def test_residue_mapping(tmp_path):
    p=tmp_path/"protein.pdb"
    p.write_text(pdb_line(1,"OD1","ASN","A",140,3,0,0,"O")+pdb_line(2,"O","HOH","A",501,2,0,0,"O","HETATM")+"END\n")
    atoms=read_receptor(p)
    assert atoms[0].residue_id=="A:ASN140"
    assert atoms[1].role=="water"

def test_exports(tmp_path):
    contacts=[ProteinContact("polar_proximity",("L:1",),("P:1",),"A:ASN140",3.0)]
    p=write_contact_csv(contacts,tmp_path/"contacts.csv")
    with p.open() as f:
        assert next(csv.DictReader(f))["residue"]=="A:ASN140"
    assert ET.fromstring(protein_contact_svg(contacts)).tag.endswith("svg")

def test_polar_and_water_screen(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    mol=Chem.MolFromSmiles("CO")
    conf=Chem.Conformer(2)
    conf.Set3D(True)
    conf.SetAtomPosition(0,Point3D(-1.4,0,0))
    conf.SetAtomPosition(1,Point3D(0,0,0))
    mol.AddConformer(conf)
    sdf=tmp_path/"ligand.sdf"
    w=Chem.SDWriter(str(sdf)); w.write(mol); w.close()
    pdb=tmp_path/"protein.pdb"
    pdb.write_text(pdb_line(1,"OD1","ASN","A",140,3,0,0,"O")+pdb_line(2,"O","HOH","A",501,1.5,2,0,"O","HETATM")+"END\n")
    contacts=screen_protein_ligand(pdb,sdf)
    assert any(c.kind=="polar_proximity" and c.residue=="A:ASN140" for c in contacts)
    assert any(c.kind=="water_bridge_candidate" and c.water_id=="P:2" for c in contacts)

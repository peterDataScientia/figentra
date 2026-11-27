"""Bond-order assignment for conventional closed-shell molecules."""
from . import Bond

def _convert(mol, ids):
    from rdkit import Chem
    kinds={Chem.BondType.SINGLE:"single",Chem.BondType.DOUBLE:"double",
           Chem.BondType.TRIPLE:"triple",Chem.BondType.AROMATIC:"aromatic"}
    out=[]
    for b in mol.GetBonds():
        if b.GetBondType() not in kinds:
            raise ValueError("Unsupported bond type: "+str(b.GetBondType()))
        out.append(Bond(ids[b.GetBeginAtomIdx()],ids[b.GetEndAtomIdx()],kinds[b.GetBondType()]))
    return out

def assign_bond_orders(atoms, *, charge, covalent_bonds=None):
    """Infer conventional bond orders from complete Cartesian geometry.

    Requires explicit hydrogens and integer total charge. Intended for
    closed-shell molecules, not transition states, metals or radicals.
    Inference may be ambiguous; prefer authoritative SDF connectivity.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from math import isfinite
    atoms=list(atoms)
    if not atoms or len({a.id for a in atoms})!=len(atoms):
        raise ValueError("Provide unique atoms")
    if not isinstance(charge,int) or isinstance(charge,bool):
        raise ValueError("Charge must be an integer")
    allowed={"H","B","C","N","O","F","Si","P","S","Cl","Br","I"}
    if any(a.element not in allowed for a in atoms):
        raise ValueError("Unsupported element for bond-order inference")
    xyz=str(len(atoms))+"\nFigentra geometry\n"+"\n".join(
        f"{a.element} {a.xyz[0]:.12g} {a.xyz[1]:.12g} {a.xyz[2]:.12g}" for a in atoms)+"\n"
    mol=Chem.MolFromXYZBlock(xyz)
    if mol is None:
        raise ValueError("Invalid XYZ geometry")
    try:
        if covalent_bonds is None:
            rdDetermineBonds.DetermineBonds(mol,charge=charge,allowChargedFragments=True)
        else:
            rw=Chem.RWMol(mol)
            ids={a.id:i for i,a in enumerate(atoms)}
            for b in covalent_bonds:
                if b.kind not in ("single","double","triple","aromatic") or b.start not in ids or b.end not in ids or b.start==b.end:
                    raise ValueError("Invalid covalent connectivity")
                rw.AddBond(ids[b.start],ids[b.end],Chem.BondType.SINGLE)
            mol=rw.GetMol()
            rdDetermineBonds.DetermineBondOrders(mol,charge=charge,allowChargedFragments=True)
        Chem.SanitizeMol(mol)
    except (RuntimeError,ValueError) as exc:
        raise ValueError("Bond-order inference failed; verify charge, hydrogens and connectivity") from exc
    if Chem.GetFormalCharge(mol)!=charge or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Inferred molecule violates charge or closed-shell assumptions")
    return _convert(mol,[a.id for a in atoms])

def read_sdf_bonds(path, *, record=0):
    """Read stored SDF connectivity and RDKit-perceived aromaticity.

    IDs are one-based atom indices matching read_sdf().
    """
    from rdkit import Chem
    supplier=Chem.SDMolSupplier(str(path),removeHs=False)
    if record<0 or record>=len(supplier) or supplier[record] is None:
        raise ValueError("Invalid SDF record")
    mol=supplier[record]
    return _convert(mol,[str(i+1) for i in range(mol.GetNumAtoms())])

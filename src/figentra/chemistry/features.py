"""Chemical feature perception from sanitized RDKit molecular graphs.

Features describe the supplied protonation state, not a predicted pH state.
Atom IDs default to one-based SDF indices; pass an explicit mapping otherwise.
"""
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

@dataclass(frozen=True)
class ChemicalFeatures:
    donor_ids: tuple
    acceptor_ids: tuple
    positive_ids: tuple
    negative_ids: tuple
    hydrophobic_ids: tuple
    aromatic_rings: tuple
    metal_ids: tuple
    coordinating_ids: tuple
    halogen_donors: tuple

    def contact_kwargs(self):
        """Arguments compatible with detect_non_covalent."""
        return {k:getattr(self,k) for k in (
            "positive_ids","negative_ids","hydrophobic_ids","aromatic_rings",
            "metal_ids","coordinating_ids")}

@lru_cache(maxsize=1)
def _factory():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/"BaseFeatures.fdef"))

def perceive_features(mol, *, atom_ids=None):
    """Perceive features using RDKit BaseFeatures and conservative graph rules.

    Positive/negative IDs use explicit formal charge only. Hydrophobic atom
    sites are neutral carbons attached only to carbon or hydrogen. Metal
    coordinating IDs are candidate N/O/S acceptors and anionic N/O/S atoms.
    Halogen donor pairs are (attached carbon ID, Cl/Br/I ID); actual sigma
    holes and coordination suitability require further validation.
    """
    from rdkit import Chem
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Provide a nonempty RDKit molecule")
    mol=Chem.Mol(mol)
    Chem.SanitizeMol(mol)
    n=mol.GetNumAtoms()
    ids=tuple(str(i+1) for i in range(n)) if atom_ids is None else tuple(atom_ids)
    if len(ids)!=n or len(set(ids))!=n or any(not isinstance(i,str) or not i for i in ids):
        raise ValueError("Provide one unique nonempty string ID per atom")
    fam={"Donor":set(),"Acceptor":set()}
    for feature in _factory().GetFeaturesForMol(mol):
        if feature.GetFamily() in fam:
            fam[feature.GetFamily()].update(feature.GetAtomIds())
    pos={a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()>0}
    neg={a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()<0}
    hyd={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6
         and a.GetFormalCharge()==0 and all(b.GetAtomicNum() in (1,6) for b in a.GetNeighbors())}
    metals={a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in
            {"Zn","Mg","Ca","Fe","Mn","Cu","Co","Ni","Cd","Na","K"}}
    coord={i for i in fam["Acceptor"]|neg if mol.GetAtomWithIdx(i).GetSymbol() in ("N","O","S")}
    rings=[]
    for ring in Chem.GetSymmSSSR(mol):
        ring=tuple(ring)
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring) and all(
            mol.GetBondBetweenAtoms(a,b).GetIsAromatic() for a,b in zip(ring,ring[1:]+ring[:1])):
            rings.append(tuple(ids[i] for i in ring))
    hal=[]
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("Cl","Br","I") and a.GetFormalCharge()==0 and a.GetDegree()==1:
            b=a.GetNeighbors()[0]
            if b.GetAtomicNum()==6:
                hal.append((ids[b.GetIdx()],ids[a.GetIdx()]))
    def mapped(indices):
        return tuple(ids[i] for i in sorted(indices))
    return ChemicalFeatures(mapped(fam["Donor"]),mapped(fam["Acceptor"]),
        mapped(pos),mapped(neg),mapped(hyd),tuple(rings),mapped(metals),
        mapped(coord),tuple(hal))

def features_from_smiles(smiles, *, atom_ids=None):
    """Perceive features for a SMILES molecular graph, without generating geometry."""
    from rdkit import Chem
    return perceive_features(Chem.MolFromSmiles(smiles),atom_ids=atom_ids)

def features_from_sdf(path, *, record=0):
    """Preserve SDF atom indices and explicit hydrogens."""
    from rdkit import Chem
    supplier=Chem.SDMolSupplier(str(path),removeHs=False)
    if record<0 or record>=len(supplier) or supplier[record] is None:
        raise ValueError("Invalid SDF record")
    return perceive_features(supplier[record])

def screen_sdf_contacts(path, *, record=0, groups=None, **thresholds):
    """Perceive features and screen one 3D SDF record without reordering atoms.

    Default groups are connected components. Cross-component contacts only.
    Returns (features, atom_pair_hydrogen_bonds, other_contacts). Hydrogen
    bonds require explicit H and are restricted to different groups.
    """
    from rdkit import Chem
    from .geometry import read_sdf
    from .bond_orders import read_sdf_bonds
    from .detection import detect_hydrogen_bonds
    from .noncovalent import detect_non_covalent
    atoms=read_sdf(path,record=record)
    bonds=read_sdf_bonds(path,record=record)
    features=features_from_sdf(path,record=record)
    if groups is None:
        mol=Chem.SDMolSupplier(str(path),removeHs=False)[record]
        groups={str(i+1):k for k,frag in enumerate(Chem.GetMolFrags(mol)) for i in frag}
    if any(a.id not in groups for a in atoms):
        raise ValueError("Group mapping must cover all atoms")
    donors=set(features.donor_ids)
    donor_h={b.end for b in bonds if b.start in donors}|{b.start for b in bonds if b.end in donors}
    hb=[b for b in detect_hydrogen_bonds(atoms,bonds,acceptor_ids=features.acceptor_ids)
        if b.start in donor_h and groups[b.start]!=groups[b.end]]
    contacts=detect_non_covalent(atoms,groups=groups,
        excluded_pairs=[(b.start,b.end) for b in bonds],**features.contact_kwargs(),**thresholds)
    return features,hb,contacts

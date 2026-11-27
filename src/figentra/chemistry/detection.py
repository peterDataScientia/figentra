"""Conservative candidate connectivity and hydrogen-bond detection.

Distances in angstroms; no bond-order, protonation or energetic inference.
For known structures prefer authoritative connectivity (e.g. SDF).
"""
from itertools import combinations
from math import isfinite
from . import Bond, bond_distance, bond_angle

# Approximate single-bond covalent radii in angstroms.
_RADII={"H":.31,"B":.84,"C":.76,"N":.71,"O":.66,"F":.57,
        "Si":1.11,"P":1.07,"S":1.05,"Cl":1.02,"Br":1.20,"I":1.39}
_CAP={"H":1,"B":3,"C":4,"N":4,"O":2,"F":1,"Si":4,
      "P":5,"S":6,"Cl":1,"Br":1,"I":1}

def _atoms(atoms):
    atoms=list(atoms)
    if len({a.id for a in atoms})!=len(atoms):
        raise ValueError("Atom IDs must be unique")
    return atoms

def detect_covalent_bonds(atoms, *, tolerance=.35, min_distance=.4):
    """Infer candidate adjacency using radii and greedy coordination limits.

    All returned edges are single-style display bonds, not assigned bond orders.
    Unsupported elements raise rather than guessing; not for metal complexes.
    """
    atoms=_atoms(atoms)
    if not isfinite(tolerance) or tolerance<0 or not isfinite(min_distance) or min_distance<=0:
        raise ValueError("Invalid distance thresholds")
    if any(a.element not in _RADII for a in atoms):
        raise ValueError("Unsupported element; supply explicit connectivity")
    candidates=[]
    for a,b in combinations(atoms,2):
        d=bond_distance(a.xyz,b.xyz)
        r=_RADII[a.element]+_RADII[b.element]
        if min_distance<=d<=r+tolerance:
            candidates.append((d/r,a.id,b.id,a,b))
    counts={a.id:0 for a in atoms}
    bonds=[]
    for _,__,___,a,b in sorted(candidates,key=lambda x:x[:3]):
        if counts[a.id]<_CAP[a.element] and counts[b.id]<_CAP[b.element]:
            bonds.append(Bond(a.id,b.id,"single"))
            counts[a.id]+=1
            counts[b.id]+=1
    return bonds

def detect_hydrogen_bonds(atoms, covalent_bonds, *, acceptor_ids,
                          max_ha=2.5,max_da=3.5,min_angle=120.):
    """Return geometric D-H...A candidates using explicit acceptor IDs.

    Explicit hydrogens and covalent connectivity are required. No periodic
    boundary handling. The caller must validate acceptor chemistry.
    """
    atoms=_atoms(atoms)
    if not all(isfinite(v) for v in (max_ha,max_da,min_angle)) or min(max_ha,max_da)<=0 or not 0<=min_angle<=180:
        raise ValueError("Invalid hydrogen bond thresholds")
    by_id={a.id:a for a in atoms}
    acceptors=set(acceptor_ids)
    if not acceptors<=by_id.keys():
        raise ValueError("Unknown acceptor ID")
    if any(by_id[i].element=="H" for i in acceptors):
        raise ValueError("Hydrogen cannot be an acceptor in this detector")
    adjacency={a.id:set() for a in atoms}
    for b in covalent_bonds:
        if b.kind not in ("single","double","triple") or b.start not in by_id or b.end not in by_id or b.start==b.end:
            raise ValueError("Expected valid covalent bonds")
        adjacency[b.start].add(b.end)
        adjacency[b.end].add(b.start)
    result=[]
    for h in atoms:
        if h.element!="H" or len(adjacency[h.id])!=1:
            continue
        d=by_id[next(iter(adjacency[h.id]))]
        if d.element not in ("N","O","F"):
            continue
        for aid in sorted(acceptors):
            a=by_id[aid]
            if aid in (h.id,d.id) or aid in adjacency[d.id]:
                continue
            ha=bond_distance(h.xyz,a.xyz)
            da=bond_distance(d.xyz,a.xyz)
            if 0<ha<=max_ha and 0<da<=max_da and bond_angle(d.xyz,h.xyz,a.xyz)>=min_angle:
                result.append(Bond(h.id,aid,"hydrogen",f"{ha:.3f} Å"))
    return result

def detect_bonds(atoms, *, acceptor_ids=None, tolerance=.35):
    """Infer candidate covalent bonds and optionally hydrogen bonds."""
    atoms=list(atoms)
    bonds=detect_covalent_bonds(atoms,tolerance=tolerance)
    if acceptor_ids is not None:
        bonds+=detect_hydrogen_bonds(atoms,bonds,acceptor_ids=acceptor_ids)
    return bonds

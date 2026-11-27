"""Geometry-based non-covalent contact candidates (angstroms/degrees).

Chemical feature membership must be supplied by the caller. No energetic,
protonation, residue assignment, or periodic boundary inference is performed.
"""
from dataclasses import dataclass
from math import isfinite, acos, degrees
from itertools import combinations
from . import bond_distance, Bond

@dataclass(frozen=True)
class Interaction:
    kind: str
    first: tuple
    second: tuple
    distance: float
    angle: float = None
    offset: float = None

_METALS={"Zn","Mg","Ca","Fe","Mn","Cu","Co","Ni","Cd","Na","K"}

def detect_non_covalent(atoms, *, positive_ids=(), negative_ids=(),
                        hydrophobic_ids=(), aromatic_rings=(), metal_ids=(),
                        coordinating_ids=(), groups=None, excluded_pairs=(),
                        salt_cutoff=4.0, hydrophobic_cutoff=4.0,
                        pi_cutoff=5.5, cation_pi_cutoff=6.0,
                        metal_cutoff=2.8, parallel_angle=30.,
                        max_pi_offset=2.5, max_cation_offset=2.5):
    """Screen explicitly supplied chemical features for contact candidates.

    groups: optional mapping atom ID -> molecule/partner ID; only cross-group
    contacts are retained. excluded_pairs: explicit atom pairs to omit, e.g.
    covalent neighbours. Ring geometry uses SVD plane fitting.
    Metal distances need element-specific adjustment for quantitative work.
    """
    import numpy as np
    atoms=list(atoms)
    by={a.id:a for a in atoms}
    if len(by)!=len(atoms):
        raise ValueError("Duplicate atom IDs")
    thresholds=(salt_cutoff,hydrophobic_cutoff,pi_cutoff,cation_pi_cutoff,
                metal_cutoff,max_pi_offset,max_cation_offset)
    if any(not isfinite(v) or v<=0 for v in thresholds) or not 0<=parallel_angle<=45:
        raise ValueError("Invalid screening thresholds")
    def ids(values):
        v=tuple(sorted(set(values)))
        if any(i not in by for i in v):
            raise ValueError("Unknown atom ID")
        return v
    pos,neg,hyd,met,coord=map(ids,(positive_ids,negative_ids,
                                  hydrophobic_ids,metal_ids,coordinating_ids))
    if set(pos)&set(neg):
        raise ValueError("Conflicting charge features")
    if any(by[i].element not in _METALS for i in met):
        raise ValueError("Unsupported metal feature")
    excluded={frozenset(p) for p in excluded_pairs}
    if groups is not None and any(i not in groups for i in by):
        raise ValueError("Group mapping must cover all atoms")
    def eligible(a,b):
        return (not set(a)&set(b) and
                not any(frozenset((i,j)) in excluded for i in a for j in b)
                and (groups is None or
                     not {groups[i] for i in a}&{groups[j] for j in b}))
    out=[]
    def pair(kind,a,b,cutoff):
        if eligible((a,),(b,)):
            d=bond_distance(by[a].xyz,by[b].xyz)
            if 0<d<=cutoff:
                out.append(Interaction(kind,(a,),(b,),d))
    for a in pos:
        for b in neg: pair("salt_bridge",a,b,salt_cutoff)
    for a,b in combinations(hyd,2):
        pair("hydrophobic",a,b,hydrophobic_cutoff)
    for a in met:
        for b in coord: pair("metal_coordination",a,b,metal_cutoff)
    rings=[]
    seen=set()
    for raw in aromatic_rings:
        ring=ids(raw)
        if len(ring)<3 or ring in seen:
            raise ValueError("Provide distinct aromatic rings with at least three atoms")
        seen.add(ring)
        xyz=np.asarray([by[i].xyz for i in ring])
        center=xyz.mean(axis=0)
        _,s,vh=np.linalg.svd(xyz-center,full_matrices=False)
        if s[1]<1e-8 or np.max(np.abs((xyz-center)@vh[-1]))>.25:
            raise ValueError("Degenerate or nonplanar aromatic ring")
        if groups is not None and len({groups[i] for i in ring})!=1:
            raise ValueError("Ring spans multiple partner groups")
        rings.append((ring,center,vh[-1]))
    for (a,ca,na),(b,cb,nb) in combinations(rings,2):
        if not eligible(a,b): continue
        delta=cb-ca
        d=float(np.linalg.norm(delta))
        if not 0<d<=pi_cutoff: continue
        angle=degrees(acos(min(1.,abs(float(na@nb)))))
        h1,h2=abs(float(delta@na)),abs(float(delta@nb))
        off=max(float(np.linalg.norm(delta-(delta@na)*na)),
                float(np.linalg.norm(delta-(delta@nb)*nb)))
        if angle<=parallel_angle and off<=max_pi_offset and min(h1,h2)>=2.5:
            out.append(Interaction("pi_parallel",a,b,d,angle,off))
        elif angle>=90-parallel_angle and max(h1,h2)>=2.5 and min(h1,h2)<=1.5:
            out.append(Interaction("pi_t_shaped",a,b,d,angle))
    for a in pos:
        for b,center,normal in rings:
            if not eligible((a,),b): continue
            delta=np.asarray(by[a].xyz)-center
            d=float(np.linalg.norm(delta))
            h=abs(float(delta@normal))
            off=float(np.linalg.norm(delta-(delta@normal)*normal))
            if 0<d<=cation_pi_cutoff and h>=2.0 and off<=max_cation_offset:
                out.append(Interaction("cation_pi",(a,),b,d,None,off))
    return out

def interaction_bonds(interactions):
    """Convert atom-pair contacts to SVG display bonds.

    Ring-centroid interactions require a dedicated centroid renderer and are
    deliberately rejected rather than misleadingly mapped to one ring atom.
    """
    out=[]
    for c in interactions:
        if len(c.first)!=1 or len(c.second)!=1:
            raise ValueError("Ring contacts require centroid rendering")
        out.append(Bond(c.first[0],c.second[0],c.kind,f"{c.distance:.3f} Å"))
    return out

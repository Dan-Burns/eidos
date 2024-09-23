from openmm import *
from openmm.app import *
import MDAnalysis as mda
import numpy as np

'''
from amber14, if you needed to specify the bond length and force constant for a disulfide
bond_force = HarmonicBondForce()
bond_length = 0.20379999999999998
force_constant = 138908.79999999996

for sulfur_pair in sulfurs_to_bond:
    bond_force.addBond(sulfur_pair[0], sulfur_pair[1], bond_length, force_constant)
system.addForce(bond_force)
'''
# TODO need simple functions to loop through topology iterables and return things based on indices, resnames, and object type
# These need to be simplified and offer a cutoff to consider for forming disulfide bonds.

def find_partner_cysteines(model):
    '''
    TODO optionally return distances to decide which ones to bond.
    '''
    topology = model.topology
    positions = model.positions
    cysteines = {}
    for res in topology.residues():
        if res.name in ["CYS","CYX"]:
            for atom in res.atoms():
                if atom.name == "SG":
                    cysteines[res.index] = atom.index
    res_indices = list(cysteines.keys())
    sg_indices = list(cysteines.values())
    sg_positions = np.array(positions)[sg_indices]
    size = len(cysteines.values())
    distances = np.zeros((size,size))
    # set self distances to infinity
    np.fill_diagonal(distances,np.inf)
    for i in range(size):
        for j in range(size):
            if i == j:
                continue
            distances[i,j] = np.linalg.norm(sg_positions[i]-sg_positions[j])._value
    mins = distances.argmin(axis=1)
    # make sure they all have unique partners
    if len(np.unique(mins)) < size:
        print("Not all cysteines have unique bonding partners. Stopping")
        return
    pairs = [tuple(sorted([res_indices[i], res_indices[j]])) for i,j in enumerate(mins)]
    pairs = list(set(pairs))
    SG_ixs = [(cysteines[i],cysteines[j]) for i,j in pairs]
       
    return pairs, SG_ixs


def get_bonded_sulfurs(topology):
    bonded_cysteines = []
    bonded_sulfurs = []
    for res in topology.residues():
        if res.name in ["CYS","CYX"]:
            for bond in res.bonds():
                if ("SG" == bond.atom1.name) and ("SG" == bond.atom2.name):
                    bonded_sulfurs.append(tuple(sorted([bond.atom1.index,bond.atom2.index])))

    for pair in bonded_sulfurs:
        cysteines = []
        for res in topology.residues():
            if res.name in ["CYS","CYX"]:
                for atom in res.atoms():
                    if atom.index in pair:
                        cysteines.append(res.index)
        bonded_cysteines.append(tuple(sorted(cysteines)))
        
    return list(set(bonded_sulfurs)), list(set(bonded_cysteines))


def get_cysteine_HG_indices(topology,cys_indices=None):
    '''
    cys_indice : list of tuples
        Specify the cysteine pairs that you want to make bonded and
        collect the HG indices from these only.
    '''
    if cys_indices is not None:
        indices = [i for tup in cys_indices for i in tup]
    atoms_to_remove = []
    for res in topology.residues():
        if res.name in ["CYS","CYX"]:
            for atom in res.atoms():
                if atom.name == "HG":
                    if (cys_indices is not None):
                        if res.index in indices:
                            atoms_to_remove.append(atom)
                    else:
                        atoms_to_remove.append(atom)
                
    return atoms_to_remove

def get_cysteine_SG_atoms(topology,cys_pairs=None):
    '''
    cys_pairs : list of tuples
        Specify the cysteine pairs that you want to make bonded and
        collect the SG indices from these only.
    '''
    SGs_to_bond = []
    for cys_pair in cys_pairs:
        atoms_to_bond = []
        for res in topology.residues():
            if (res.name in ["CYS","CYX"]) and (res.index in cys_pair):
                for atom in res.atoms():
                    if atom.name == "SG":
                        atoms_to_bond.append(atom)
        SGs_to_bond.append(atoms_to_bond)
                
    return SGs_to_bond

def get_pairs_to_bond(bonded, pairs):
    return [pair for pair in pairs if pair not in bonded]

def make_disulfide_bonds(model, SGs_to_bond=None):
    '''
    Parameters
    ----------
    model : openmm.app.modeller
    '''
    topology = model.topology
    cys_pairs, SG_pairs = find_partner_cysteines(model)
    bonded_sulfurs, bonded_cysteines = get_bonded_sulfurs(topology)
    sulfur_pairs_to_bond = get_pairs_to_bond(bonded_sulfurs, SG_pairs)
    cysteine_pairs_to_bond = get_pairs_to_bond(bonded_cysteines, cys_pairs)
    HGs_to_remove = get_cysteine_HG_indices(topology,cysteine_pairs_to_bond)
    model.delete(HGs_to_remove)
    topology=model.topology
    if SGs_to_bond is not None:
        for pair in SGs_to_bond:
        topology.addBond(*pair)
    else:
        SGs_to_bond = get_cysteine_SG_atoms(topology, cysteine_pairs_to_bond)
        for pair in SGs_to_bond:
            topology.addBond(*pair, 'Single')
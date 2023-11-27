## trjconvert scripts 
## path to gromacs bin
import parmed as pmd
from parmed.gromacs import GromacsTopologyFile
import copy

def parmed_underscore_topology(gromacs_processed_top, atom_indices, output_top):
    '''
    Add underscores to atom types of selected atoms.
    This is useful if using the plumed_scaled_topologies script 
    for hremd system modification.
    With this, you still need to open the new topology file and delete the 
    underscores from the beginning of the file [atomtypes]
    or else plumed will look for atoms with 2 underscores to apply lambda to.
    '''
    top = GromacsTopologyFile(gromacs_processed_top)

    for atom in top.view[atom_indices].atoms:
        atom.type = f"{atom.type}_"
        if atom.atom_type is not pmd.UnassignedAtomType:
            atom.atom_type = copy.deepcopy(atom.atom_type)
            atom.atom_type.name = f"{atom.atom_type.name}_"


    top.save(output_top)

def get_plumed_distance_restraints(atom_dictionary, universe):
    '''
    Convenience function to print out pairs of atom indices in plumed
    .dat file format and return a dictionary of interatomic distances
    to use in specifying distance restraints.
    
    Parameters
    ----------
    atom_dictionary : dictionary
        Key and values are tuples of ("resname", "residue id", "atom name")
        identifying atom pairs to apply distance restraints to
    
    universe : mda.Universe
        The universe containing the starting structure with distances to base 
        restraitns off of.

    Returns
    -------
    Prints out the format to copy and paste into a plumed file to specify the atom pairs
    for plumed to track the distances of.

    Returns a dictionary of distances for the identifiers in the printed results.

    Example
    -------
    pairs = {
        ("1PR", "10", "O'L"): ("Ser", "225" ,"HG"),
        ("1PR", "10", "O8B",1): ("Ser", "225" ,"H")
        }

    distances = get_plumed_distance_restraints(pairs, u)
    <out> 0: DISTANCE ATOMS=7232,3451
          1: DISTANCE ATOMS=7226,3442
    print(distances)
    <out> {0: 1.6858517,
          1: 1.8405167}
    
    Plumed file:
    UNITS LENGTH=A
    0: DISTANCE ATOMS=7232,3451
    1: DISTANCE ATOMS=7226,3442

    UWALL1: UPPER_WALLS ...
        ARG=0,1
        AT=4,4    # added extra space
        KAPPA=100,100
        ...

          
    
    '''
    u = universe
    distances = {}
    for i, (key, value) in enumerate(pairs.items()):
        a = u.select_atoms(f'resname {key[0]} and resnum {key[1]} and name {key[2]}').atoms[0]
        a_id = a.id
        a_pos = a.position
        b = u.select_atoms(f'resname {value[0].upper()} and resnum {value[1]} and name {value[2]}').atoms[0]
        b_id = b.id
        b_pos = b.position
        dist = np.linalg.norm(a_pos - b_pos)
        distances[i]=dist
        print(f'{i}: DISTANCE ATOMS={a_id},{b_id}')
    return distances
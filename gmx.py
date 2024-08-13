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

################### In progress - make function to automatically select the group to use in processing
################### You might want to center on a ligand but output the system minus water
################## alternatively, you could run make_ndx and save the group menu and cancel
################# Then construct run it again and make a group from the components of interest
################## Then you'll know ahead of time how many menus a given command will return and 
################### use the info dictionary to provide the centering, clustering, and output groups.....
import subprocess
import os
def get_gromacs_groups(gmx_folder, index_file=None, verbose=False):
    '''
    Uses gmx make_ndx to retrieve the indices of the different system components.
    You can use these to automatically make selections in gmx commands that 
    have interactive prompts.

    gmx_folder : str
        Path to the folder containing the gromacs files. The folder and files
        are automatically produced in omm.OMMSetup.save().
    
    index_file : bool
        If True, the make_ndx command will include this file and the returned 
        groups will include the ids of any additional groups in the index file
        that are not automatically available from the .gro file.

    Returns
    -------
    List of tuples (group_index: group_name)
    # TODO the info dictionary should have the something to connect the omm
    # chains to the gromacs groups so that you can automatically produce a
    # gromacs index file that has a group containing all of the system components
    # you want in a reduced post-processed trajectory e.g. protein and ligand only.

    
    '''
    try:
        # get the .gro file
        files = os.listdir(gmx_folder)
        gro_file = [f'{gmx_folder}/{file}' for file in files if file.endswith('gro')][0]
        
        # Start the gmx command to get the available groups
        if index_file == True:
            index_file = [f'{gmx_folder}/{file}' for file in files if file.endswith('ndx')][0]
            command = ['gmx', 'make_ndx', '-f', gro_file, '-n', index_file]
        else:
            command = ['gmx', 'make_ndx', '-f', gro_file]
        process = subprocess.Popen(
            command,  
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True  # Make sure to handle text mode
        )
        
        try:
            # Wait for the process to complete or timeout
            stdout, stderr = process.communicate(timeout=.1) # how small of a timeout?
        except subprocess.TimeoutExpired:
            # If the process times out, kill it
            process.kill()
            stdout, stderr = process.communicate()
        
        # Print the raw output for debugging
        if verbose:
            print("Raw output:\n", stdout)
        
        # Parse the group list from the output
        groups = []
        for line in stdout.split('\n'):  
            line = line.lstrip().split()
            if len(line)>1:
                if line[0].isdigit():
                    group_number = line[0]
                    group_name = line[1]
                    groups.append((group_number, group_name))
        
        return groups

    except Exception as e:
        print(f"An error occurred: {e}")
        print("stdout : ", stdout)
        print("stderr : ", stderr)
        return []

################# use system/simulation_info.json ################
def run_trjconv(group_dict):
    # Use pty to handle the pseudo-terminal for interactivity
    master, slave = pty.openpty()
    
    # Start the trjconv process
    process = subprocess.Popen(
        ['gmx', 'trjconv', '-s', 'input.tpr', '-f', 'input.xtc', '-o', 'output.xtc'],
        stdin=slave,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True  # Make sure to handle text mode
    )
    
    # Close the slave end of the pty
    os.close(slave)
    
    # Read the output to get the available groups
    group_list = ""
    while True:
        line = os.read(master, 1024).decode('utf-8')
        group_list += line
        if "Select a group:" in line:
            break
    
    # Print the captured group list for debugging
    print("Available groups:\n", group_list)
    
    # Parse the available groups and match with the dictionary
    group_selection = None
    for line in group_list.split('\n'):
        for key in group_dict:
            if key in line:
                group_selection = group_dict[key]
                break
    
    if group_selection is not None:
        # Send the selected group
        os.write(master, f"{group_selection}\n".encode('utf-8'))
    else:
        raise ValueError("No matching group found in the output.")
    
    # Read the remaining output
    stdout, stderr = process.communicate()
    
    # Print the output and error (if any)
    print(stdout)
    if stderr:
        print("Error:", stderr)

# Example group dictionary
group_dict = {
    "Protein": "1",
    "Water": "2"
}

run_trjconv(group_dict)

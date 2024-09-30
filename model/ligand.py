import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Bond
import mendeleev as mv
import parmed


class MolFromPDB():
    '''
    UNTESTED 

    Take in a pdb and make edits so a resulting .mol file
    will have charges and bond types specified.  This is needed
    if you are going to use openmmforcefields to parameterize a
    small molecule.
    TODO allow for smiles input and manipulation as well

    pdb : str
        Path to pdb file.  The pdb should have the correct protonation state.

    
    '''

    # parmed df makes it easy to find the atom name from rdkits GetIDX

    def __init__(self, pdb, sanitize=True, removeHs=False ):
        self.mol = Chem.MolFromPDBFile(pdb, sanitize=sanitize, removeHs=removeHs)
        # get 2d representation
        #Chem.AllChem.Compute2DCoords(self.mol)


    def mol_with_atom_index(self):
        #TODO display large and readable - looks incomprehensible as is
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return self.mol
    
    def set_bond_type(self, idx1, idx2, bondtype):
        '''
        Change existing bond type to bondtype

        Example
        -------
        double_bonds = [("O2B", "C2B"), ("O9", "C9"), ("O'1", "P7'"), ("O'L", "P7B")]

        for pair in double_bonds:
            indices = []
            for at in mol.mol.GetAtoms():
                name = lig_df.iloc[at.GetIdx()]['name'] # lig_df is parmed structure to_dataframe()
                if name in pair:
                    indices.append(at.GetIdx())
            print(f'Setting {indices} to double bonded.')
            mol.set_bond_type(indices[0], indices[1], 'double')
        '''
        bondtypes = {'single':Chem.rdchem.BondType.SINGLE,
                     'double':Chem.rdchem.BondType.DOUBLE,
                     'triple':Chem.rdchem.BondType.TRIPLE}


        for bond in self.mol.GetBonds():
            if bond.GetEndAtomIdx() in {idx1, idx2} and bond.GetBeginAtomIdx() in {idx1, idx2}:
                bond.SetBondType(bondtypes[bondtype])

    def set_charge(self, idxs, charge):
        '''
        
        Example
        -------
        negative_charges = ["O8'", "O9'", "O8B", "O9B"]

        indices = []
        for atom in negative_charges:
            for at in mol.mol.GetAtoms():
                name = lig_df.iloc[at.GetIdx()]['name'] # lig_df is parmed structure to_dataframe()
                if name == atom:
                    indices.append(at.GetIdx())
                    print(f'Setting {at.GetIdx()} to charge of -1.')
        mol.set_charge(indices, -1)
        '''

        for at in self.mol.GetAtoms():
            if at.GetIdx() in idxs: #and mv.element(at.GetAtomicNum()).name == 'Oxygen':
                at.SetFormalCharge(charge)

    def show_mol_file(self):
        print(Chem.MolToMolBlock(self.mol))

    def save_mol_file(self, output):
        if output.split('.')[-1] != 'mol':
            output = f'{output}.mol'
        Chem.MolToMolFile(self.mol,output)

def make_PET_polymer_smiles(n_units, end='ethylene_glycol'):
    '''
    For working with PET simulations. Construct PET polymers of n_units.

    n_units : int
        The number of Terephthalate monomers
    
    Returns
    -------
    Smiles string for PET polymer of n_units.

    Example
    --------
    smiles = make_PET_polymer_smiles(4)
    molecule = Chem.MolFromSmiles(smiles)
    molecule = AllChem.AddHs(molecule)
    # remove the hydrogens on the carboxy termini if desired
    remove_carboxy_hydrogens(molecule)

    '''

    monomer = "OC(=O)C1=CC=C(C=C1)C(=O)O"
    l_term = "OC(=O)C1=CC=C(C=C1)C(=O)"
    r_term = "C(=O)C1=CC=C(C=C1)C(O)=O"
    internal = "C(=O)C1=CC=C(C=C1)C(=O)"
    linker = "OCCO"
    if n_units == 1:
        if end == 'ethylene_glycol':
            return l_term+linker
        else:
            return monomer
    elif n_units == 2:
        if end == 'ethylene_glycol':
            return l_term+linker+internal+linker
        else:
            return l_term+linker+r_term
    else:
        polymer = l_term
        for unit in range(n_units-2):
            addition = linker+internal
            polymer+=addition
        if end == 'ethylene_glycol':
            return polymer+=linker+internal+linker
        else:
            return polymer+=linker+r_term

def find_carboxy_hydrogens(molecule):
    carboxy_termini_hydrogens = []
    for atom in molecule.GetAtoms():
        # Carbons bonded to 2 oxygens
        if atom.GetSymbol() == "C" and (len([neighbor.GetIdx() for neighbor in atom.GetNeighbors()
                                            if neighbor.GetSymbol() == 'O'])==2): 
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    carboxy_termini_hydrogens.extend([Oneighbor.GetIdx() for Oneighbor in neighbor.GetNeighbors() 
                                                     if Oneighbor.GetSymbol() == 'H'])
    
    return carboxy_termini_hydrogens
    
def remove_carboxy_hydrogens(PET_molecule):
    '''
    Provide a protonated RDKit Mol object of a PET molecule and 
    remove the hydrogens on the termini so they are negatively charged.

    PET_molecule : RDKit Mol object of PET

    Returns
    -------
    RDKit Mol
    '''

    ixs = find_carboxy_hydrogens(PET_molecule)
    

    edit_mol = Chem.RWMol(PET_molecule)
    for ix in ixs:
        edit_mol.RemoveAtom(ix)

    clean_mol = edit_mol.GetMol()
    for atom in clean_mol.GetAtoms():
        if (atom.GetAtomicNum() == 8) and (len([bond.GetIdx() for bond in atom.GetBonds()]))==1:
            if "SINGLE" in [bond.GetBondType().name for bond in atom.GetBonds()]:
                atom.SetFormalCharge(-1)
    Chem.SanitizeMol(clean_mol)
    clean_mol = AllChem.AddHs(clean_mol)
    return clean_mol

############################## Functions to explore setting up grid points for RESP fitting ###############################
def center_molecule(u):
    '''
    Parameters
    ----------
    u : mda.Universe
        Universe containing only the atoms that will be translated

    Returns
    -------
    mda.Universe
        The universe where the atoms have been translated to be centered at the origin
    '''
    centroid = np.mean(u.atoms.positions, axis=0)
    centered_coords = u.atoms.positions - centroid
    u.atoms.positions = centered_coords
    return u

def get_min_max_dims(u):
    '''
    Parameters
    ----------
    u : mda.Universe
        Universe containing only the atoms that will be translated

    Returns
    -------
    np.array, np.array
    The first array contains the minimum coordinates, the second array the maximum
    '''
    return np.min(u.atoms.positions,axis=0), np.max(u.atoms.positions,axis=0)
    
def generate_rectangular_mesh_points(min_dims, max_dims, spacing=0.2):
    '''
    Generate a rectangular prism of mesh points
    
    '''
    x = np.arange(min_dims[0], max_dims[0], spacing)
    y = np.arange(min_dims[1], max_dims[1], spacing)
    z = np.arange(min_dims[2], max_dims[2], spacing)

    X,Y,Z = np.meshgrid(x, y, z)

    points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]) # 3 X N points (x points, y points, z points)
    return points.T


def fibonacci_sphere(atom_name, atom_coords, table, n_points=20, extension=1.4):
    """
    Generates exactly n_points uniformly distributed on a sphere using the Fibonacci lattice method.
    Can use to generate efficient number of grid points for RESP fitting. Each call returns an array of points 
    arranged in a sphere around the given point so after all the points have been generated, points less than the vdW 
    distance from any atom in the molecule must be removed. 

    Parameters:
    atom_name : str
        The chemical symbol for the atom

    atom_coords : np.ndarray
        The coordinates for the atom
    
    table : RDKit.Chem.GetPeriodicTable()
        Instantiated RDKit PeriodicTable object

    n_points : int
        Number of points to generate.

    extension : float
        The distance beyond the van der Waals radius to generate the points.
    
    Returns:
    np.ndarray: Array of shape (n_points, 3) containing the generated points.
    """
    radius = table.GetRvdw(atom_name) + extension
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians

    for i in range(n_points):
        y = 1 - (i / float(n_points - 1)) * 2  # y goes from 1 to -1 (latitude)
        r = np.sqrt(1 - y * y)  # Radius at y
        theta = phi * i  # Increment based on the golden angle

        x = np.cos(theta) * r
        z = np.sin(theta) * r
        points.append([x * radius, y * radius, z * radius])

    return np.array(points)

def generate_random_sphere_points(atom_name, atom_coords, table, n_points=20, extension=1.4):
    """
    Generates n_points uniformly distributed on a sphere of given radius.
    
    Parameters:
    n_points (int): Number of points to generate.
    radius (float): The radius of the sphere.
    
    Returns:
    np.ndarray: Array of shape (n_points, 3) containing the generated points.
    """
    radius = table.GetRvdw(atom_name) + extension
    points = []
    for _ in range(n_points):
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        x = radius * np.sin(phi) * np.cos(theta)
        y = radius * np.sin(phi) * np.sin(theta)
        z = radius * np.cos(phi)
        points.append([atom_coords[0] + x, atom_coords[1] + y, atom_coords[2] + z])
    
    return np.array(points)

def remove_points_within_vdw(u, points, table, extension):
    """
    Removes points that fall within the van der Waals radius of any atom in the molecule.


    Parameters:
    u : mda.Universe
    points (np.ndarray): Array of shape (n_points, 3) containing the coordinates of the points.
    table (RDKit PeriodicTable): RDKit PeriodicTable instance to retrieve van der Waals radii.

    Returns:
    np.ndarray: Array of points that are outside the van der Waals radii of all atoms.
    """

    # Initialize a boolean mask for points (True means the point is kept)
    mask = np.ones(len(points), dtype=bool)

    # Loop through each atom in the molecule
    for atom in u.atoms:
        atom_index = atom.ix
        atom_name = atom.name
        atom_pos = atom.position
        radius = table.GetRvdw(atom_name) + extension  # Get the van der Waals radius of the atom
        
        # Calculate the distance from the atom to each point
        distances = np.linalg.norm(points - atom_pos, axis=1)
        
        # Mark points that are within the vdW radius for removal
        mask &= (distances > radius)
    
    # Return only the points that are outside the vdW radius of all atoms
    return points[mask]

def generate_grid_points(u, table, n_points=50, extension=1.4, random=False):

    '''
    Universe should already be centered


    random : bool
        If True, generate points with np.random.uniform. 
        Otherwise distribute points evenly with on a fibonacci sphere.

    Returns
    -------
    Dictionary of atom indices and np.arrays of points in a sphere around each atom

    Visualize the grid points with nglview:
    view = nglview.show_mdanalysis(u)
    view
    for point in fixed_points:
        view.shape.add_sphere(point,[0,0,1],0.2)
    '''
    grid_points = {}
    
   
    # Loop through each atom in the molecule
    for atom in u.atoms:
        atom_index = atom.ix
        atom_name = atom.name
        atom_pos = atom.position
        
        # Get the van der Waals radius of the atom and extend it
        vdw_radius = table.GetRvdw(atom_name)  # Default to 1.5 if not found
        radius = vdw_radius + extension
        
        # Generate the grid points for the atom
        if random==True:
            sphere_points = generate_random_sphere_points(atom_name, atom_pos, table=table, n_points=n_points, extension=extension)
        else:
            sphere_points = fibonacci_sphere(atom_name, atom_pos, table=table, n_points=n_points, extension=extension)
        
        # Translate the points to be centered around the atom's position
        # translated_points = sphere_points + np.array([atom_pos.x, atom_pos.y, atom_pos.z])
        
        # Store the points associated with this atom in the dictionary
        grid_points[atom_index] = sphere_points

    # stack all the oints together
    points = np.vstack([array for array in grid_points.values()])
    # remove any points that fall within the vdW radius 
    fixed_points = remove_points_within_vdw(u, points, table, extension)
    
    return fixed_points

    ######################################### End grid point section #################################################


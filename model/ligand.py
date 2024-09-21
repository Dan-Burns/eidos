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

def make_PET_polymer(n_units):
    '''
    For working with PET simulations. Construct PET polymers of n_units.

    n_units : int
        The number of Terephthalate monomers
    
    Returns
    -------
    Smiles string for PET polymer of n_units.
    '''
    monomer = "OC(=O)C1=CC=C(C=C1)C(=O)O"
    l_term = "OC(=O)C1=CC=C(C=C1)C(=O)"
    r_term = "C(=O)C1=CC=C(C=C1)C(O)=O"
    internal = "C(=O)C1=CC=C(C=C1)C(=O)"
    linker = "OCCO"
    if n_units == 1:
        return monomer
    elif n_units == 2:
        return l_term+linker+r_term
    else:
        polymer = l_term
        for unit in range(n_units-2):
            addition = linker+internal
            polymer+=addition
        polymer+=linker+r_term
        return polymer
    
def remove_carboxy_hydrogens(PET_molecule):
    '''
    Provide a protonated RDKit Mol object of a PET molecule and 
    remove the hydrogens on the termini so they are negatively charged.

    PET_molecule : RDKit Mol object of PET

    Returns
    -------
    RDKit Mol
    '''

    ixs = []
    for atom in PET_molecule.GetAtoms():
        if atom.GetAtomicNum() == 1 and ("O" in [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]):
            ixs.append(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                for bond in neighbor.GetBonds():
                    print(bond.GetIdx())
    ixs.sort(reverse=True)

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

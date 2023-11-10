import rdkit
from rdkit import Chem
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

    def __init__(self, pdb, removeHs=False ):
        self.mol = Chem.MolFromPDBFile(pdb,removeHs=removeHs)


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


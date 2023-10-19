import rdkit
from rdkit import Chem
from rdkit.Chem import Bond
import mendeleev as mv

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
    def __init__(self, pdb, removeHs=False ):
        self.mol = Chem.MolFromPDBFile(pdb,removeHs=removeHs)


    def mol_with_atom_index(self):
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return self.mol
    
    def set_bond_type(self, idx1, idx2, bondtype):
        
        bondtypes = {'single':Chem.rdchem.BondType.SINGLE,
                     'double':Chem.rdchem.BondType.DOUBLE,
                     'triple':Chem.rdchem.BondType.TRIPLE}


        for bond in self.mol.GetBonds():
            if bond.GetEndAtomIdx() in {idx1, idx2} or bond.GetBeginAtomIdx() in {idx1, idx2}:
                bond.SetBondType(bondtypes[bondtype])

    def set_charge(self, idxs, charge):

        for at in self.mol.GetAtoms():
            if at.GetIdx() in idxs: #and mv.element(at.GetAtomicNum()).name == 'Oxygen':
                at.SetFormalCharge(charge)

    def show_mol_file(self):
        print(Chem.MolToMolBlock(self.mol))

    def save_mol_file(self, output):
        if output.split('.')[-1] is not 'mol':
            output = f'{output}.mol'
        Chem.MolToMolFile(self.mol,output)
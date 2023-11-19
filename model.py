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



from openff.toolkit.topology import *
import openmmtools as omt
from openmmtools.integrators import *
from openforcefields import *

import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import numpy as np


def minimize_morph(structure, reference, output):
    '''
    force a structure into the backbone conformation of another reference structure.
    Useful for homologs where you might have an experimental structure and then a modeled structure with a different sequence.
    If you think the modelled structure is more likely to take the conformation of the experimental structure,
    you can make the morph the modelled structure into the backbone conformation of the experimental structure.
    
    '''

    ### NOT TESTED, Just ported it over from a notebook 

    structure = mda.Universe(structure)

    reference = mda.Universe(structure)

    # identical number of residues?
    n_residues = int(structure.residues.n_residues)
    reference.residues.n_residues == structure.residues.n_residues

    # make a ndarray to populate with the distances
    distance_array = np.zeros((n_residues, n_residues))

    # How far do you want to search for neighboring C alphas?
    cutoff = 10
    # Go through all the atoms and find the CAs and the neighboring CAs
    for atom in u_aws.atoms:
        if atom.name == 'CA':
            neighbors = structure.select_atoms(f'around {cutoff} index {atom.ix} and name CA')
            for neighbor_atom in neighbors.atoms:
                
                if neighbor_atom.name == 'CA' and neighbor_atom.ix != atom.ix:
                
                    distance_array[atom.residue.ix][neighbor_atom.residue.ix] = np.linalg.norm(
                                                                atom.position-neighbor_atom.position)
    
    # get the equilibrium distances for the restraint
    r0_distances = distance_array[distance_array>4]

    # get the residues to apply potential to
    residue_indices = np.argwhere(distance_array>4)

    ###############################################################################
    ## Simulation Set Up - can use the minimize sidechains function above
    # use implicit solvent since we're just minimizing and forcing into conformation
    pdb = PDBFile(structure)
    # have to adjust this for implicit solvent 
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
   
    modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer, ionicStrength=0.1*molar )
    temperature = 303.15*kelvin
    integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
    pressure = 1 * bar
    system = forcefield.createSystem(modeller.topology,nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
            constraints=HBonds)
    system.addForce(MonteCarloBarostat(pressure, temperature))
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    #################################################################################

    ########### Restraints ################
    # map the residue indices to their CA atom indices
    CA_res_atom_indices = {}
    for atom in list(modeller.topology.atoms()):
        if atom.name == 'CA':
            CA_res_atom_indices[atom.residue.index] = atom.index

    # restraint_force[0].setBondParameter

    restraint_weight = 10000 * openmm.unit.kilocalories_per_mole / openmm.unit.angstrom ** 2
    restraint_force = CustomBondForce("0.5*K*(r-r0)^2")
    # Add the restraint weight as a global parameter in kcal/mol/nm^2
    restraint_force.addGlobalParameter("K", restraint_weight)
    # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
    restraint_force.addPerBondParameter("r0")
    custom_forces ={}
    for i, pair in enumerate(residue_indices):
        
        restraint_force.addBond(CA_res_atom_indices[pair[0]],CA_res_atom_indices[pair[1]], [r0_distances[i]])
    custom_forces['distance_restraints']=system.addForce(restraint_force)
    # necessary?
    simulation.context.reinitialize(preserveState=True)


    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True)

    with open(output, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.positions, f)


from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import app
from openmm import unit



            


#### Don't Touch This file - Seq2Ensemble imports this

def minimize_sidechains(output, pdb_file, temperature=300.00):

    '''
    Energy minimizes a protein structure with backbone restraints.
    Used after modelling sidechain atoms back onto a coarse grained structure.

    output: str
        file name to write the minimized structure to.

    pdb_file: str
        path to pdb structure file to be minimized
    
    temperature: float or int
        temperature in K for minimization
    '''
    # load the pdb, probably could have skipped Modeller() and just used pdb.topology and pdb.positions
    pdb = PDBFile(pdb_file)
    modeller = Modeller(pdb.topology, pdb.positions)

    temperature = temperature*kelvin

    forcefield = ForceField('amber14-all.xml','implicit/obc2.xml')
    system = forcefield.createSystem(modeller.topology,nonbondedMethod=app.NoCutoff,
                                    constraints=app.HBonds, 
                                    
    )
   
    # using implicit solvent, no cutoff is used
    #removed this argument for openmm 8.0 -- implicitSolvent=app.OBC2,
    '''
    openmm 7 version
    forcefield = ForceField('amber14-all.xml')
    system = forcefield.createSystem(modeller.topology,nonbondedMethod=app.NoCutoff,
                                    constraints=app.HBonds, implicitSolvent=app.OBC2,
                                    implicitSolventSaltConc=0.1*moles/liter,
    )
    '''
    ## CREATE THE SIMULATION
    integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
    # No pressue with implicit solvent
    #pressure = 1 * bar
    #system.addForce(MonteCarloBarostat(pressure, temperature))
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    ### GET RESTRAINED ATOMS
    topology = simulation.topology
    # atoms to restrain
    restraint_indices = []
    # get the backbone CA, N, C
    for atom in topology.atoms():
        if atom.name == 'CA' or \
            atom.name == 'N' or \
            atom.name == 'C':
            restraint_indices.append(atom.index)



    # Add position restraints to heavy atoms to allow water to relax around protein
    # Create the restraint object, force, and add the particles to it
    positions = simulation.context.getState(getPositions=True).getPositions()
    reference_coordinates = positions.in_units_of(unit.nanometer)
    restraint_weight = 500 * unit.kilocalories_per_mole / unit.angstrom ** 2
    restraint_force = CustomExternalForce('K*periodicdistance(x, y, z, x0, y0, z0)^2')
    # Add the restraint weight as a global parameter in kcal/mol/nm^2
    restraint_force.addGlobalParameter("K", restraint_weight)
    # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    for index in range(0, len(positions)):
        if index in restraint_indices:
            xyz = reference_coordinates[index].in_units_of(unit.nanometers) / unit.nanometers
            restraint_force.addParticle(index, xyz)
    custom_forces ={}
    custom_forces['positional_restraints'] = system.addForce(restraint_force)


    ##################################### END Create Position Restraints #################

    ####### Energy Minimize #####################
    simulation.context.reinitialize()
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    with open(output, 'w') as f:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
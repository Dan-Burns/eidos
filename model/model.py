from openff.toolkit.topology import *
import openmmtools as omt
from openmmtools.integrators import *
from openforcefields import *

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import app
from openmm import unit


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

############################### Sali Lab Modeller Functions #############################################################3

from Bio import pairwise2, SeqIO
import Bio
from modeller import *
import re
from modeller.automodel import *
# Salilab Modeller functions for missing loops

def make_aln_file(pdb_path, output):
    '''
    pdb_path : string 
        Path to pdb file.

    output : string
        File path ending in ".seq" extension.
    '''
    # https://salilab.org/modeller/wiki/Missing_residues

    code = re.split(r'\.|/',pdb_path)[-2]
    e = Environ()
    m = Model(e, file=pdb_path)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=output)

def write_modeller_alignment(seqa, seqb, output, seqa_header = ['seqa','None'], seqb_header = ['seqb','None']):
    with open(output, 'w') as f:
        f.write(f">P1;{seqa_header[0]}\n{seqa_header[1]}\n{seqa}*\n")
        f.write(f">P1;{seqb_header[0]}\n{seqb_header[1]}\n{seqb}*\n")

def get_modeller_alignment(pdb_path, seq_file, output):
    '''
    seq_file : string   
        Path to file from make_aln_file.
    '''
    code = re.split(r'\.|/',pdb_path)[-2]

    for record in SeqIO.parse(pdb_path, "pdb-seqres"):
        print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        print(record.dbxrefs)
    
        # get full seq including missing residues
        # https://biopython.org/wiki/SeqIO
        # https://salilab.org/modeller/8v2/manual/node176.html (pir format)
    for data in SeqIO.parse(seq_file, "pir"):
        print(data)
        if type(data) == Bio.SeqRecord.SeqRecord:
            pdb_seq = data
    # https://biopython.org/docs/1.75/api/Bio.pairwise2.html
    ## USING localmx here
    alignments = pairwise2.align.localmx(pdb_seq.seq, record.seq, 1, 0)
    seqa, seqb = alignments[0].seqA, alignments[0].seqB
    # PIR format 
    # First line is "">P1;pdb"
    # fields for second line
    # 1 specification whether 3d structure is available and method used to obtain "structure" is enough
    # 2 PDB code or filename (if no path is provided modeller searches IOData.atom_files_directory)
    # 3-6 Residue and chain identifiers first(fields 3-4) last (fields 5-6)
    # 7 optional protein name
    # 8 optional source
    # 9 optional resolution
    # 10 optional R-factor
    write_modeller_alignment(seqa, seqb, output, seqa_header = [code,f'structureX:{code}:FIRST:A:LAST:A::::'], 
                        seqb_header = [f'{code}_fill', 'sequence:::::::::'])


def fill_loops(structure_dir, ali_file, pdb_code):
    # ali file needs 
    # >P1;1abc
    # structureX:>P1;6yhr_fill:FIRST:A:LAST:A::::
    # >P1;1abc_fill
    # for the two seqs assuming pdb code is 1abc

    log.verbose()
    env = Environ()

    # directories for input atom files
    env.io.atom_files_directory = ['.', structure_dir]

    a = LoopModel(env, alnfile = ali_file,
                knowns = pdb_code, sequence = f'{pdb_code}_fill')
    # might have to go in to the ali_file and adjust the chain ID
    # or remake the aln file
    # this is sensitive to the way the pdb is setup - multiple chains complicates things
    a.starting_model= 1
    a.ending_model  = 1

    a.loop.starting_model = 1
    a.loop.ending_model   = 2
    a.loop.md_level       = refine.fast

    a.make()

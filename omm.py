import openmm
from openmm import *
from openmm.app import *
from openmm.vec3 import Vec3
from openmm.unit import *
import MDAnalysis as mda
from MDAnalysis.transformations.rotate import rotateby
import nglview as nv
import numpy as np
from pdbfixer import PDBFixer
import os
import subprocess
import shutil
import parmed
import json


def omm_to_mda(topology,positions):
    '''
    Provide an OpenMM topology and positions and return an MDAnalysis Universe.
    
    Parameters
    ----------
    
    Returns
    -------
    mda.Universe
    '''
    top = topology
    positions
    topology = {}
    topology['n_atoms'] = len(list(top.atoms()))
    topology['n_residues'] = len(list(top.residues()))
    topology['n_segments'] = len(list(top.chains()))
    topology['atom_resindex'] = [atom.residue.index for atom in list(top.atoms())] #resids for each atom
    topology['residue_segindex'] = [res.chain.index for res in list(top.residues())] #segids for each residue
    u = mda.Universe.empty(**topology,trajectory=True)
    u.add_TopologyAttr('name', [atom.name for atom in top.atoms()])
    u.add_TopologyAttr('type', [atom.element.symbol for atom in top.atoms()])
    u.add_TopologyAttr('resname', [res.name for res in top.residues()])
    u.add_TopologyAttr('resid', [res.index for res in top.residues()])
    u.add_TopologyAttr('segid', [chain.id for chain in top.chains()])
    u.add_TopologyAttr('mass', [atom.element.mass.value_in_unit(dalton) for atom in top.atoms()])
    #TODO add mass
    u.atoms.positions = np.asarray([np.asarray(position._value) for position in positions.in_units_of(angstroms)])
    bonds = []
    for bond in list(top.bonds()):
        bonds.append((bond.atom1.index,bond.atom2.index))
    u.add_TopologyAttr('bonds', bonds)
    u.add_TopologyAttr("chainID")
    for seg in u.segments:
        sel = u.select_atoms(f'segid {seg.segid}')
        sel.atoms.chainIDs = f'{seg.segid}'
 
    return u

def top_pos_from_sim(simulation):
    state = simulation.context.getState(getPositions=True)
    return simulation.topology, state.getPositions()

def make_vec3(array):
    '''
    take a 3 component array (a vector) and turn into a Vec3
    '''
    return Vec3(*array)

def omm_to_pdb(omm_object, output):
    '''
    Provide an OpenMM object that contains topology and positions (Modeller,PDBFile, or Simulation)
    and write it to a pdb file.
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    
    if type(omm_object) == openmm.app.modeller.Modeller or\
        type(omm_object) == openmm.app.pdbfile.PDBFile:
        top, pos = omm_object.topology, omm_object.positions   

    elif type(omm_object) == openmm.app.simulation.Simulation:
        top, pos = top_pos_from_sim(omm_object) 
    with open(output, 'w') as f:
        PDBFile.writeFile(top, pos, file=f)



class PositionModifier():
    '''
    provide an openmm model and the interactively modify the positions of the constituent molecules.
    Translations and rotations are visualized with nglview.
    The prinicipal axes are shown.
    Best to use just the solutes and move them around before solvating.
    '''
    def __init__(self, model):
        self.model = model
        self.top, self.original_pos = model.topology, model.positions
        self.u = omm_to_mda(self.top, self.original_pos)

    def translate(self, selection, direction, distance):
        '''
        distance: float or int
            distance to translate the selection in agstroms
        '''
        u = self.u.copy()
        
        directions = {'x':[1,0,0], 'y':[0,1,0], 'z':[0,0,1]}
        if type(direction) == str and direction in directions.keys():
            direction = np.array(directions[direction])
            
        else:
            try: 
                np.reshape(direction,3)
            except:
                print('Expecting a string (x, y, or z) or a vector.')
                
        index = np.where(direction == 1)[0][0]
        
        if type(selection) ==  str:
            selection = u.select_atoms(selection)
        translation_array = np.zeros((selection.atoms.positions).shape)
        translation_array[:,index] = distance
        selection.atoms.positions += translation_array
        self.u = u

    def rotate(self, selection, direction, angle):
        '''
        
        '''
        u = self.u.copy()
        sel = u.select_atoms(selection)
        cog = sel.center_of_geometry()
        ts = sel.ts
        new_coords = rotateby(angle, direction, point=cog, ag=sel)(ts)
        new_coords.positions
        sel.atoms.positions = new_coords.positions
        self.u = u

    def view(self, axes=True):
        #https://nglviewer.org/ngl/api/manual/selection-language.html
        #https://nglviewer.org/ngl/api/manual/index.html
        # axes just stretches out a single axis as you move the pieces apart
        view = nv.show_mdanalysis(self.u)
        if axes:
            view.add_axes()
        return view

    def to_model(self):
        # changes the original model
        data = list(map(make_vec3,list(self.u.atoms.positions)))
        # quantity is imported from unit
        self.model.positions = quantity.Quantity(data,unit=angstroms) 
        return self.model

    def get_closest_distance(sel1,sel2):
        sel1 = self.u.select_atoms(sel1)
        sel2 = self.u.select_atoms(sel2)
        
        return distance_array(sel1.atoms, sel2.atoms,).min()


def get_force_index(force_name, system):
    for i, force in enumerate(system.getForces()):
        if force.getName() == force_name:
            return i

def remove_force_by_name(force_name, system):
    force_index = get_force_index(force_name, system)
    system.removeForce(force_index)

def slow_heat(simulation, start_temp=1, end_temp=293, 
              nsteps=10e5):
    '''
    todo : accept Quantity class from openmm unit 
    todo : display clocktime
    Heat simulation in 1 kelvin intervals from start_temp to end_temp

    Parameters
    ----------
    simulation : openmm.app.Simulation

    start_temp : Temperature to begin heating in kelvin.

    end_temp : Temperature to stop heating in kelvin.

    nsteps : Number of steps over which the heating will occur. 
            Number of steps per kelvin is nsteps / (end_temp - start_temp)

    Returns
    -------
    Simulation with integrator temperature set to end_temp.
    '''
    steps_per_interval = int(np.floor(nsteps/(end_temp - start_temp))) 
    for temp in range(start_temp, int(end_temp+1)):
        #TODO make this print interval smart
        if temp % 25 == 0:
            print(f'Current Temperature : {temp}')
        simulation.integrator.setTemperature(temp)
        simulation.step(steps_per_interval)

#TODO add pdb2pqr

def fix_pdb(input_pdb, output_pdb, pH=7.0, keep_water=True, replace_nonstandard_resis=True):
    '''
    PDBFixer convenience function 
    
    '''
    # https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    if replace_nonstandard_resis:
        fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keep_water)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))

def get_sim_info(simulation):
    '''
    Returns a dictionary of basic simulation information.
    '''
    top = simulation.topology
    sys = simulation.system
    integrator = simulation.integrator
    forces = simulation.getForces()
    info = {}
    disulfide_bond_list = []
    # get disulfide bonds
    for bond in top.bonds():
        if bond.atom1.name == 'SG' and bond.atom2.name == 'SG':
            disulfide_bond_list.append(bond)
    info['disulfide_bonds'] = disulfide_bond_list
    info['num_atoms'] = top.getNumAtoms()
    info['num_chains'] = top.getNumChains()
    info['num_residues'] = top.getNumResidues()
    info['residue_names'] = set([res.name for res in top.residues()])
    info['current_step'] = simulation.currentStep
    info['integrator'] = integrator.__class__
    info['forces'] = forces
    for force in forces:
        if force.getName() == 'NonbondedForce':
            info['nonbonded_cutoff'] = force.getCutoffDistance()
    info['periodic_boundaries'] = sys.usesPeriodicBoundaryConditions()
    if sys.usesPeriodicBoundaryConditions() == True:
        info['periodic_box_vectors'] = sys.getDefaultPeriodicBoxVectors()
    info['unit_cell_dimensions'] = top.getUnitCellDimensions()
    info['step_size'] = integrator.getStepSize()
    info['friction'] = integrator.getFriction()
    info['temperature'] = integrator.getTemperature()
    info['seed'] = integrator.getRandomNumberSeed()
    
    return info


class OMMSetup:
    '''
    ################# IN PROGRESS #############################
    Class to piece together an openmm simulation object

    todo: write json file with simulation parameters that can be used to 
    load the integrator and so on for restarts
    '''

    def __init__(self, structures,
                 structures_to_parameterize=None,
                 nonbonded_cutoff=1*nanometer,
                 integrator_type=LangevinMiddleIntegrator,
                 forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
                 temperature=None,
                 pressure=1*bar,
                 box_shape='cube',
                 padding=0.6*nanometer,
                 ):
        self.structures = structures
        self.structures_to_parameterize = structures_to_parameterize
        self.nonbonded_cutoff = nonbonded_cutoff
        self.integrator_type = integrator_type
        self.forcefields = forcefields
        self.temperature = temperature
        self.pressure = pressure
        self.box_shape = box_shape
        self.padding = padding
        # add padding or box vectors, ions, concentration, water model

    '''
    structures : dict
        dict of keys of user supplied names and values of paths to prepared PDB 
        files for each component of the system.
        Example
        -------
        structures = ['lysozyme':'./253L.pdb']

    structures_to_parameterize : dictionary
        If structures contains small molecules or chemicals not available the specified
        forcefields, then supply a dictionary with names matching those provided in 
        "structures" and a value corresponding to the path to a .mol file.

    '''
    def model(self):
        # modeler components
        pdb_file = self.structures[0]
        pdb = PDBFile(pdb_file)
        modeller = Modeller(pdb.topology, pdb.positions)
        if len(self.structures) > 1:
            for structure in self.structures[1:]:
                pdb_file = structure
                pdb = PDBFile(pdb_file)
                modeller.add(pdb.topology, pdb.positions)
        self.modeller = modeller
            

    def parameterize(self):
        self.forcefield = ForceField(*self.forcefields)
        if self.structures_to_parameterize is not None:
            molecules = []
            for key, val in self.structures_to_parameterize.items():
                molecule = molecule.from_file(val)
                molecule.assign_partial_charges('am1bcc')
                molecules.append(molecule)
            smir = SMIRNOFFTemplateGenerator(molecules=molecules) # register a list?
            self.forcefield.registerTemplateGenerator(smir.generator)
        self.modeller.addSolvent(self.forcefield, padding=self.padding,
                            ionicStrength=0.1*molar, model='tip3p',
                            boxShape=self.box_shape)
    
    def make_system(self):
        # create system object
        system = self.forcefield.createSystem(self.modeller.topology, 
                                         nonbondedMethod=PME,
                                         nonbondedCutoff=self.nonbonded_cutoff, 
                                        constraints=HBonds)
        # Add pressure control
        system.addForce(MonteCarloBarostat(self.pressure, self.temperature))
        self.system=system

    def make_simulation(self):
        integrator = self.integrator_type(self.temperature, 1/picosecond, 2*femtoseconds) # add options to init
        # create simulation object
        simulation = Simulation(self.modeller.topology, self.system, integrator)
        simulation.context.setPositions(self.modeller.positions)
        simulation.minimizeEnergy()
        self.simulation = simulation


    def save(self, output, name='system'):
        '''
        save gromacs files and openmm system files

        output : str
            Path to output. Directory will be created if none exists.

        name : str
            Optional name to prefix to your saved files.
        '''
        # save the system and minimized structure
        topology, positions = top_pos_from_sim(self.simulation)
        os.makedirs(f'{output}',exist_ok=True)
        with open(f'{output}/{name}_system.xml', 'w') as outfile:
            outfile.write(XmlSerializer.serialize(self.system))
        with open(f'{output}/{name}_minimized.pdb', 'w') as f:
            PDBFile.writeFile(topology, positions, f)
        
        # create a .json file with the simulation's details to be used for 
        sim_info = get_sim_info(self.simulation)
        sim_info['forcefields'] = self.forcefields
        sim_info['box_shape'] = self.box_shape
        sim_info['padding'] = self.padding
        #sim_info['PME'] = True/False
        
        with open('sim_info.json','w') as h:
            json.dump(sim_info, h)
  
        ##### Save a gromacs topology for future trjconv use - Use a no-constraints version of system to avoid parmed error
        parmed_system = self.forcefield.createSystem(self.simulation.topology, nonbondedMethod=PME,nonbondedCutoff=1*nanometer, rigidWater=False)
        pmd_structure = parmed.openmm.load_topology(self.simulation.topology, system=parmed_system, xyz=positions)

        os.makedirs(f'{output}/gmx/',exist_ok=True)
        pmd_structure.save(f"{output}/gmx/{name}_gmx.top", overwrite=True)
        pmd_structure.save(f"{output}/gmx/{name}_gmx.gro", overwrite=True)

        
        # write an energy minimization .mdp file to use with gmx grompp
        from .files.text import em_mdp
        with open(f'{output}/gmx/em.mdp','w') as w:
            w.write(em_mdp)
        # Check to see if gmx is available
        command_path = shutil.which('gmx')
        # write a .tpr file that can be used for things like trjconv
        if command_path is not None:
            command = [
            'gmx', 'grompp',
            '-f', f'{output}/gmx/em.mdp',
            '-c', f"{output}/gmx/{name}_gmx.gro",
            '-p', f"{output}/gmx/{name}_gmx.top",
            '-o', f"{output}/gmx/{name}_em_gmx.tpr"
            ]

            # Execute the command
            result = subprocess.run(command, capture_output=True, text=True)

            # Print the output and error (if any)
            print("Output:\n", result.stdout)
            print("Error:\n", result.stderr)
    

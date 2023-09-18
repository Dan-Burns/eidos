import openmm
from openmm import System
from openmm.app import PDBFile, Modeller, Simulation
from openmm.vec3 import Vec3
from openmm.unit import *
import MDAnalysis as mda
from MDAnalysis.transformations.rotate import rotateby
import nglview as nglview
import numpy as np



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
    #TODO add mass
    u.atoms.positions = np.asarray([np.asarray(position._value) for position in positions.in_units_of(angstroms)])
    bonds = []
    for bond in list(top.bonds()):
        bonds.append((bond.atom1.index,bond.atom2.index))
    u.add_TopologyAttr('bonds', bonds)
    u.add_TopologyAttr("chainID")
    for seg in u.segments:
        sel = u.select_atoms(f'segid {seg}')
        sel.atoms.chainIDs = f'{seg}'
 
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
        
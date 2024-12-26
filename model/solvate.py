from openmm import *
from openmm.app import *
from scipy.spatial import cKDTree

# Not Tested / pseudocode atm
# tools to add custom solvent/ particles to a simulation box

pdb = PDBFiles(input_structure)
modeller = Modeller(pdb.topology, pdb.positions)

solute_positions = pdb.positions.value_in_unit(unit.nanometers)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.topology.setUnitCellDimensions((x,y,z))

custom_particle_spacing = spacing * unit.nanometers
custom_particle_sigma = .2

solute_tree = cKDTree(solute_positions)
# Generate candidate positions on a grid
x_range = np.arange(x_dim_min, x_dim_max, custom_particle_spacing(unit.nanometers))
y_range = np.arange(y_dim_min, y_dim_max, custom_particle_spacing(unit.nanometers))
z_range = np.arange(z_dim_min, z_dim_max, custom_particle_spacing(unit.nanometers))
candidate_positions = np.array(np.meshgrid(x_range, y_range, z_range)).T.reshape(-1,3)

custom_particle_positions = []
for candidate in candidate_positions:
    if not solute_tree.query_ball_point(candidate, custom_particle_sigma):
        custom_particle_positions.append(Vec3(*candidate) * unit.nanometers)

# add the particle elements to the modeller
elements = [Element.getBySymbol(el)] * len(custom_particle_positions)

# create a system with solute
system = forcefield.createSystem(modeller.topology...)
custom_chain = modeller.topology.addChain()
for position in custom_particle_positions:
    res = modeller.topology.addResidue("res", custom_chain)
    modeller.topology.addAtom("atom", element, res)
    modeller.positions.append(position)
    system.addParticle(element_mass)

for force in system.getForces():
    if force.getName() == "NonbondedForce": # assuming just adding noble gas particles - more complicated if bonded terms
        nonbonded_force = force
for i in range(len(custom_particle_positions)):
    nonbonded_force.addParticle(custom_particle_charge, custom_particle_sigma, custom_particle_epsilon)

#proceed to simulation instance
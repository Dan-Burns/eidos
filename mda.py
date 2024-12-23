import MDAnalysis as mda
from MDAnalysis.lib.util import convert_aa_code as conv_aa
from MDAnalysis.analysis import align, rms


# TODO full script or class for post processing that will spit out all the 
## vital stats, rmsd, rmsf, rgyr, sasa, pca, contacts, etc. with options to include 
## specific sub-selections like domains or active sites and visualization scripts.

def write_trajectory(u, selection='protein', output_file='saved_traj.xtc', start=0, stop=-1, step=1):
    '''
    
    '''
    sel = u.select_atoms(selection)

    sel.write(output_file, frames=u.trajectory[start:stop:step])
    # old
    '''
    with mda.Writer(output_file, sel.n_atoms) as w:
        for frame in u.trajectory[start:stop:step]:
            w.write(sel)
    '''
def get_residue_names(u):
    return list(set(u.residues.resnames))

def align_traj(u, ref=None, selection='name CA', output_file=None):
    '''
    Align a trajectory to a reference based on selection.

    u : mda.Universe
    
    ref : mda.Universe
        The reference universe. If None, align to first frame of trajectory

    selection : str
        MDAnalysis selection of atoms to base the alignment on.

    output_file : str
        Path to output trajectory if saving to file.

    Returns
    -------
    Modifies the universe in place. If output_file is provided, writes the aligned
    trajectory to output_file.
    '''
    if ref == None:
        ref = u
    if output_file is not None:
        align.AlignTraj(u, ref, select=selection, filename=output_file).run()
    else:
        align.AlignTraj(u, ref, select=selection, in_memory=True).run()

def get_rmsd(u, ref=None, selection='name CA'):
    '''
    Only basic rmsd call. 

    #TODO add more functionality

    Returns
    -------
    Nx3 np.ndarray [[frame, time (ps), RMSD (A)]]
    '''

    output = rms.RMSD(u, reference=ref, select=selection,).run() # parallelizable

    return output.rmsd

def get_sequence(u,print_seq=True,by_chain=False):
    '''
    Return the amino acid sequence of a universe.
    '''

    resis = list(map(conv_aa,[res.resname for res in u.residues]))
    # TODO make sure this just outputs the sequence of the polymers and not the solvent in the event this is a simulation system.
    resis_by_chain = {seg.segid:list(map(conv_aa,[res.resname for res in seg.residues])) for seg in u.segments}
    if print_seq:
        if by_chain:
            for chain, reslist in resis_by_chain.items():
                print(f'{chain} :',"".join(reslist))
        else:
            print("".join(resis))
    else:
        if by_chain:
            return resis_by_chain
        else:
            return resis
        

        
def process_traj(structure, trajectory,
                 sel='not (resname HOH or resname SOL or resname Na or resname Cl or resname WAT or resname Na+ or resname Cl-)',
                 reduced_traj_output = 'reduced_traj.xtc', reduced_structure='reduced_structure.pdb',
                 align=True, align_sel='name CA',
                 output_file='aligned_traj.xtc'):
    '''
    Take an input trajectory and return the trajectory aligned to a reference and minus a selection like solvent.

    Parameters
    ----------
    structure : string (path to structure)
                Pdb file matching the trajectory's topology.

    trajectory : string (path to trajectory) 
                molecular dynamics trajectory compatible with MDAnalysis

    remove_sel : string
                MDAnalysis selection of components to remove from the processed trajectory.
                If you want to keep everything and just align, set remove_sel=None.
    
    reduced_traj_output : string
                        Path to save the trajectory of selected atoms.

    reduced_structure : string
                    Path to save the structure containing the selected atoms.

    align : boolean
            
    align_sel : Selection to use as the alignment reference.  If "None" trajectory will be aligned
                to starting structure's C-alpha atoms.

    output_file : string
                    Path to output processed trajectory.  If "None", assign the in-memory
                    output universe to a variable.
    '''
    if sel is not None:
        u = mda.Universe(structure, trajectory)
        sel = u.select_atoms(sel)
        sel.write(reduced_structure)

        write_trajectory(u, selection=sel, output_file=reduced_traj_output, start=0, stop=-1, step=1)

        reduced_u = mda.Universe(reduced_structure, reduced_traj_output)
        ref = mda.Universe(reduced_structure)

        if output_file is not None and align==True:
            align.AlignTraj(reduced_u,  
                        ref,  
                        select=align_sel,  
                        filename=output_file,  
                    ).run()
        elif output_file is None and align==True:
            align.AlignTraj(reduced_u,  
                        ref,  
                        select=align_sel,  
                        in_memory=True,  
                    ).run()
            return reduced_u
    else:
        u = mda.Universe(structure, trajectory)
        ref = mda.Universe(structure)
        if output_file is not None and align==True:
            align.AlignTraj(u,  
                        ref,  
                        select=align_sel,  
                        filename=output_file,  
                    ).run()
        else:
            align.AlignTraj(u,  
                        ref,  
                        select=align_sel,  
                        in_memory=True,  
                    ).run()
            # TODO check
            # no need to return the u if it is updated in place
            

# todo Class for simulation run/continuation/processing
# use a json file for everything that updates with new current run # etc.
# include checks to make sure naming doesn't overlap - don't 
# require any specifications of run numbers etc
# handle pre heating start.

import numpy as np

def compute_angle_ACF(angles, max_tau, step=1, sample_origin='independent', type='angle'):
    '''
    Compute the autocorrelation function avg(cos[theta(t)-theta(tau+t)])tau
    see: van der Spoel, D. & Berendsen, H. J. Molecular dynamics simulations of Leu-enkephalin in water and DMSO. 
            Biophys. J. 72, 2032–2041 (1997)

    Parameters
    ----------
    angles : np.array
        M x N where M is the number of trajectory frames the angles are taken from
        and N is the number of dihedral angles.
        Angles are expected to be in degrees.
    max_tau : int
        The maximum time lag (in frames) to evaluate the function to.
    step : int
        The tau interval.  If step = 2 then the correlation function is evaluated at tau = 0, 2, 4, 6, 8 .....
    sample_origin : string 
        options are
        'all' : function is evaluated for origin time (t) of every frame (or row) in angles
        'independent' : function is evaluted for frames (or rows) equal to 0, tau, ...., ((M/tau)-1)tau
        #TODO intermediate origins
    type : string
        options are
        'angle' : inputs are expected to be angles (degrees). cos(angles[sample_origin] - angles[sample_origin+tau]) is computed.
        'vector' : inputs are expected to be vectors (taken from the cross product of two bond orientation vectors)
                    dot(vector[sample_origin], vector[sample_origin+tau]) is computed
    '''
    
    angles = np.radians(angles)
    # get the number of trajectory frames
    n_rows = len(angles)
    # create an empty array that will have a number of rows
    # equal to the number of time lag samples
    # and a number of columns equal to the number of angles 
    ACF = np.zeros((int(max_tau/step), angles.shape[1]))
    # first row contains cos of 0 = 1 (angle minus itself i.e. no lag)
    ACF[0,:] = np.ones(angles.shape[1])
    
    if step > 1:
        start=step
    else:
        start=1
    # for each time lag starting from "step" frame lags to the max number of lags in "step" size
    for i, tau in enumerate(range(start,max_tau+1,step)):
        # make an array that will contain a number of rows equal to the length of trajectory minus the time lag
        # because the final data point will be tau frames back from the final frame
        # eg 10 total rows and a tau of 2, there will be 8 rows to average because the final data point is cos(row 8 - row 10)
        if sample_origin == 'independent':
            n_samples = int((n_rows/tau)-1)
            diffs = np.zeros((n_samples, angles.shape[1]))
            coef = tau
        else:
            n_samples = int(n_rows-tau)
            diffs = np.zeros((n_samples, angles.shape[1]))
            coef = 1
        # loop over n_samples origins 
        for record, t in enumerate(range(n_samples)):
            sample_origin = t*coef
            # get the cos of the difference between angle at time t and t+tau
            if type == 'angle':
                diffs[record,:] = np.cos(angles[sample_origin] - angles[sample_origin+tau])
            elif type == 'vector':
                diffs[record,:] = np.dot(angles[sample_origin],angles[sample_origin+tau]) 
        # after they've all been recorded, get the means for each angle and
        # record them on the row of ACF corresponding to this iterations tau value
        ACF[i,:] = diffs.mean(axis=0)
    return ACF

def get_average_angle(angles):
    '''
                      sum_i_from_1_to_N sin(a[i])
    a = arctangent ---------------------------
                      sum_i_from_1_to_N cos(a[i])
    '''
    avg_angle = np.arctan(np.sin(angles.sum(axis=0)/np.cos(angles).sum(axis=0)))

    return avg_angle

def scale_to_unit(vector):
    
    return vector/np.linalg.norm(vector)

def get_bond_vector_cross_products(ag,atoms=('H','N','CA')):
    '''
    ag : mda.Univsere AtomGroup containing the atoms from the residues of interest.
    '''
    cross_products = []
    residues = ag.residues
    # get the row indices corresponding to the atoms
    i = np.where(residues.atoms.names == atoms[0])[0]
    j = np.where(residues.atoms.names == atoms[1])[0]
    k = np.where(residues.atoms.names == atoms[2])[0]
    # get the position vectors for the atoms
    i,j,k = residues.atoms.positions[[i,j,k]].squeeze()
    # cross product
    n = scale_to_unit(np.cross(i-j, j-k))
    cross_products.append(n)
    return cross_products
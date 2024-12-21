import MDAnalysis as mda
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
import numpy as np
import pandas as pd

# https://userguide.mdanalysis.org/stable/examples/analysis/custom_trajectory_analysis.html
class RadiusOfGyration(AnalysisBase):  # subclass AnalysisBase

    def __init__(self, atomgroup, verbose=True):
        """
        Set up the initial analysis parameters.
        """
        # must first run AnalysisBase.__init__ and pass the trajectory
        trajectory = atomgroup.universe.trajectory
        super(RadiusOfGyration2, self).__init__(trajectory,
                                               verbose=verbose)
        # set atomgroup as a property for access in other methods
        self.atomgroup = atomgroup
        # we can calculate masses now because they do not depend
        # on the trajectory frame.
        self.masses = self.atomgroup.masses
        self.total_mass = np.sum(self.masses)

    def _prepare(self):
        """
        Create array of zeroes as a placeholder for results.
        This is run before we begin looping over the trajectory.
        """
        # This must go here, instead of __init__, because
        # it depends on the number of frames specified in run().
        self.results = np.zeros((self.n_frames, 6))
        # We put in 6 columns: 1 for the frame index,
        # 1 for the time, 4 for the radii of gyration

    def _single_frame(self):
        """
        This function is called for every frame that we choose
        in run().
        """
        # call our earlier function
        rogs = radgyr(self.atomgroup, self.masses,
                      total_mass=self.total_mass)
        # save it into self.results
        self.results[self._frame_index, 2:] = rogs
        # the current timestep of the trajectory is self._ts
        self.results[self._frame_index, 0] = self._ts.frame
        # the actual trajectory is at self._trajectory
        self.results[self._frame_index, 1] = self._trajectory.time

    def _conclude(self):
        """
        Finish up by calculating an average and transforming our
        results into a DataFrame.
        """
        # by now self.result is fully populated
        self.average = np.mean(self.results[:, 2:], axis=0)
        columns = ['Frame', 'Time (ps)', 'Radius of Gyration',
                   'Radius of Gyration (x-axis)',
                   'Radius of Gyration (y-axis)',
                   'Radius of Gyration (z-axis)',]
        self.df = pd.DataFrame(self.results, columns=columns)
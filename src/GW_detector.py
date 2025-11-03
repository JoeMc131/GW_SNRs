
import numpy as np
from scipy.interpolate import interp1d
import sys

# list the detectors
detectors = [
    'CE',
    'ET',
    'LISA'
]

#paths to the sensitivity curves
paths = [
    '../data/ce_strain/cosmic_explorer_40km_lf_strain.txt',
    '../data/ET-0000A-18_ETDSensitivityCurveTxtFile.txt',
    '../data/LISA_sensitivity_curves.txt'
]

class GW_detector:
    def __init__(self, detector):

        if detector not in detectors:
            print('\nError: Detector must be one of ["LISA", "ET", "CE"]')
            sys.exit(0)

        for i in range(len(detectors)):

            # read the sensitivity curve data_columns into the object
            if detector == detectors[i]:

                self.data = np.loadtxt(paths[i])
                self.freq_max = max(self.data[:, 0])
                self.freq_min = min(self.data[:, 0])

                self.name = detectors[i]

                self.freq = self.data[:,0]

                # the strain data_columns is in the fourth column for the ET data_columns,
                # for the rest its the first column
                if detector == 'ET':
                    self.strain = self.data[:,3]
                else:
                    self.strain = self.data[:,1]

                self.characteristic_strain = self.sensitivity_curves()

                break


    def sensitivity_curves(self):

        interpolated_line = interp1d(self.freq, self.strain)

        return interpolated_line

import GW_detector as gw
import Signal as si

import numpy as np
import matplotlib.pyplot as plt

import sys

solar_mass = 1.989e30

def plot_merger(detector, merger):
    """
    Plots the strain of the detector and the merger signal

    :param detector:
    :param merger:
    :return: nothing
    """

    # calculate the SNR
    snr = merger.calc_SNR(detector)

    if snr == 'n/a':
        snr = 0

    # plot figure
    plt.figure()
    plt.xlabel('Frequency [Hz]', size = 15)
    plt.ylabel('Strain [Hz$^{-1/2}$]', size = 15)
    plt.plot(detector.freq, detector.strain, label = f'{detector.name} Characteristic Strain')
    plt.plot(merger.f, merger.strain(), label = 'Signal')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.legend()
    plt.title(f'SNR = {snr:.1e}')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    # record the inputted arguments
    m1 = float(sys.argv[1])
    m2 = float(sys.argv[2])
    z = float(sys.argv[3])
    detector_arg = sys.argv[4]

    # set up the detector
    detector = gw.GW_detector(detector_arg)

    # set up the merger
    merger = si.Signal(m1 = m1*solar_mass, m2 = m2*solar_mass, z = z)


    #calculate the SNR
    snr = merger.calc_SNR(detector)

    if snr == 'n/a':
        snr = 0


    ############### OUTPUT #####################

    print('\n################# SUMMARY ###################'
          '\n'
          f'\nMASS 1:      {m1:.1e}'
          f'\nMASS 2:      {m2:.1e}'
          f'\nREDSHIFT:    {z}'
          f'\nSNR:         {snr:.1e}'
          f'\n'
          f'\n#############################################')

    # plot the merger
    plot_merger(detector, merger)

"""
Code adapted from Robson, Cornish and Liu 2018 SNR calc prescription

"""

import numpy as np
from astropy.cosmology import WMAP7 as cosmo

a_k = [2.9740e-1,5.9411e-1,5.0801e-1,8.4845e-1]
b_k = [4.4810e-2,8.9794e-2,7.7515e-2,1.2848e-1]
c_k = [9.5560e-2,1.9111e-1,2.2369e-2,2.7299e-1]


G = 6.674e-11
c = 299792458

year = 365.25*24*60*60


def chirp_mass(m1, m2):
    M = m1 + m2
    return ((m1*m2)**(3/5))/(M**(1/5))


class Signal:
    def __init__(self, m1, m2, z, T_obs = 4*year):
        self.Mass_1 = m1
        self.Mass_2 = m2
        self.redshift = z

        self.T_obs = T_obs

        self.chirp_mass = chirp_mass(self.Mass_1, self.Mass_2)
        self.eta = (self.Mass_1*self.Mass_2)/(self.Mass_1 + self.Mass_2)**(2)

        self.f0 = self.f_k(0)
        self.f1 = self.f_k(1)
        self.f2 = self.f_k(2)
        self.f3 = self.f_k(3)

        a = (G*self.chirp_mass/c**3)

        self.f_start = (1/(8*np.pi*a))*(5*a/(T_obs))**(3/8)
        self.f_end = self.f3

        self.f = np.logspace(np.log10(self.f_start),
                             np.log10(self.f_end),
                             1000)


    def f_k(self, k):
        M = self.Mass_1 + self.Mass_2
        fk = (a_k[k]*self.eta**2 + b_k[k]*self.eta + c_k[k])
        fk/=(np.pi * (G*M)/(c**3))
        return fk

    def L(self, f):
        ans = (1/(2*np.pi))*(self.f2/((f - self.f1)**2 + (self.f2**2)/4))
        return ans

    def w(self):
        ans = (np.pi*self.f2)/2
        ans *= (self.f0/self.f1)**(2/3)
        return ans

    def A(self, f):
        A_0 = np.sqrt(5/24)
        A_0 *= (G*self.chirp_mass/c**3)**(5/6) * self.f0**(-7/6)
        A_0 /= ((np.pi**(2/3) * (cosmo.luminosity_distance(self.redshift).to('m').value)/c))

        mask1 = (f < self.f0)
        mask2 = (f >= self.f0) & (f < self.f1)
        mask3 = (f >= self.f1) & (f < self.f3)

        A_array = np.zeros(len(f))

        A_array[mask1] = A_0 * (f[mask1]/self.f0)**(-7/6)
        A_array[mask2] = A_0 * (f[mask2]/self.f0)**(-2/3)
        A_array[mask3] = A_0 * self.w() * self.L(f[mask3])

        return A_array

    def calc_SNR(self, detector):

        detector_f = detector.freq
        gw_f = self.f


        # cut of the signal frequency domain that does not lie in the detector frequency domain
        if self.f_start < min(detector_f):
            mask = gw_f >= min(detector_f)
            gw_f = gw_f[mask]
        if self.f_end > max(detector_f):
            mask = gw_f <= max(detector_f)
            gw_f = gw_f[mask]


        if(len(gw_f)==0):
            return 0

        # # set integration limit
        # mask = (detector_f >= self.f_start) & (detector_f <= self.f_end)
        # detector_f = detector_f[mask]
        else:
            f = np.logspace(np.log10(min(gw_f)), np.log10(max(gw_f)), 100000)

            N = len(f)

            d_ln_f = np.log(f[1:]) - np.log(f[:N-1])

            integrand_1 = (f[1:]*self.A(f)[1:])/detector.sensitivity_curves()(f[1:])
            integrand_2 = (f[:N-1]*self.A(f)[:N-1])/detector.sensitivity_curves()(f[:N-1])

            snr = np.sum(0.5*(integrand_1 + integrand_2)*d_ln_f)

            return np.sqrt(16/5)*snr

    def strain(self):
        return np.sqrt(self.f*self.A(self.f)**2)






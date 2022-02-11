"""
file name: isrf.
author: Sacha Gavino
date: Nov 2021
language: PYTHON 3.8
_____________________________________________________________________________________________________
short description:  create ISRF as a function of wavelengths.
_____________________________________________________________________________________________________
"""
import numpy as np

class InterstellarRadFields:
    def __init__(self, cut=2.e-1, d78=True, vdb82=True):
        self.cut=cut 
        self.d78=d78 
        self.vdb82=vdb82


    def draine78(self, lam):
        '''
        Desc: Draine (1978) fit from Sterberg & Dalgarno (1995) between 912 and 2000 Angstrom
        Args: lam (microns)
        return: 2D array (wavelength [nm], ISRF [erg.cm-2.s-1.hz-1.sr-1])
        '''
        #lam_d78 = self.lam[self.lam <= 2.e-1]*1e3 #take the right range and convert to nm.
        d78 = (3.2028e13*lam**(-3) - 5.1542e15*lam**(-4) + 2.0546e17*lam**(-5) ) / (4*np.pi) #divide by 4pi to get sr-1
        d78 *= (2.998e-12)/(2.998e17) # 1photon = 3e-12 ergs and 1nm = 2.998e17 Hz
        d78 = np.stack([lam, d78])
        return d78
        

    def van_black82(self, lam):
        '''
        Desc: Extension of van Dishoek & Black (1982) at longer wavelengths (> 2000 Angstrom)
        Args: lam (microns)
        '''
        #lam_vdb82 = self.lam[self.lam > 2.e-1]*1e3
        vdb82 = 3.67e4*lam**(0.7) / (4*np.pi)
        vdb82 *= (2.998e-12)/(2.998e17)
        vdb82 = np.stack([lam, vdb82])
        return vdb82

    def create_isrf(self, lam):
        lam = lam*1e3
        if self.d78 is True and self.vdb82 is True:
            lam_d78 = lam[lam <= self.cut*1e3]
            lam_vdb82 = lam[lam > self.cut*1e3]
            isrf = np.concatenate((self.draine78(lam_d78), self.van_black82(lam_vdb82)), axis=1)

        elif self.d78 is True and self.vdb82 is False:
            isrf = self.draine78(lam) 

        elif self.d78 is False and self.vdb82 is True:
            isrf = self.van_black82(lam)

        return isrf




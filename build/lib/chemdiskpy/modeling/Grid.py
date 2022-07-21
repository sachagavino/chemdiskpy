"""
_____________________________________________________________________________________________________________
file name: Grid
@author: P.Sheehan. Adapted by S. Gavino for chemistry codes.
last update: Aug 2021
language: PYTHON 3.8
short description:  class Grid for young stellar objects modeling. 
_____________________________________________________________________________________________________________
"""

import numpy as np

class Grid:
    def __init__(self):
        self.density = []
        self.dustdensity = []
        self.gasdensity_chem = []
        self.dustdensity_chem = []
        self.dustdensity_single_chem = []
        self.hg_chem = []
        self.tgas_chem = []
        self.temperature = []
        self.localfield = []
        self.avz = []
        self.stars = []
        self.isrf = []
        self.dust = []
        self.accretionheating = []

    def add_star(self, star):
        self.stars.append(star)

    def add_isrf(self, isrf):
        self.isrf.append(isrf)

    def add_temperature(self, temperature):
        self.temperature.append(temperature)

    def add_localfield(self, localfield):
        self.localfield.append(localfield)

    def add_density(self, density):
        self.density.append(density)

    def add_dustdensity(self, density):
        self.dustdensity.append(density)

    def add_dustdensity_chem(self, density):
        self.dustdensity_chem.append(density)

    def add_dustdensity_single_chem(self, density):
        self.dustdensity_single_chem.append(density)

    def add_gasdensity_chem(self, density):
        self.gasdensity_chem.append(density)

    def add_gastemperature_chem(self, gas_temperature):
        self.tgas_chem.append(gas_temperature)

    def add_hg_chem(self, hg):
        self.hg_chem.append(hg)

    def add_avz(self, av_z):
        self.avz.append(av_z)

    def add_dust(self, dust):
        self.dust.append(dust)

    def add_accretionheating(self, q_visc):
        self.accretionheating.append(q_visc)
        
    def set_cartesian_grid(self, xmin, xmax, nx):
        #w1, w2, w3 provide grid with coordinates using the center of each cell.
        self.coordsystem = "cartesian"

        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(xmin, xmax, nx)
        z = np.linspace(xmin, xmax, nx)

        w1 = 0.5*(x[0:x.size-1] + x[1:x.size])
        w2 = 0.5*(y[0:y.size-1] + y[1:y.size])
        w3 = 0.5*(z[0:z.size-1] + z[1:z.size])

        return np.stack((x, y, z)), np.stack((w1, w2, w3))


    def set_spherical_grid(self, w1, w2, w3):
        self.coordsystem = "spherical"

        self.r = 0.5*(w1[0:w1.size-1] + w1[1:w1.size])
        self.theta = 0.5*(w2[0:w2.size-1] + w2[1:w2.size])
        self.phi = 0.5*(w3[0:w3.size-1] + w3[1:w3.size])

        self.w1 = w1
        self.w2 = w2
        self.w3 = w3

    def set_nautilus_grid(self, r, h_lim=4, nz_chem=64):
        #hg = self.disk.scaleheight(np.array(r))
        pts = np.arange(0, nz_chem, 1)
        #zchem = np.ones((len(rchem), nb_points))

        #hh, ptpt = np.meshgrid(hchem, pts)
        z = (1. - (2.*pts/(2.*nz_chem - 1.)))*h_lim#*Hg
        self.rchem = np.array(r)
        self.zchem = z
        self.nz_chem = nz_chem

    def set_wavelength_grid(self, lmin, lmax, nlam, log=False): #microns
        if log:
            self.lam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.lam = np.linspace(lmin, lmax, nlam)

    def set_mcmonowavelength_grid(self, lmin, lmax, nlam, log=False):
        if log:
            self.monolam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.monolam = np.linspace(lmin, lmax, nlam)

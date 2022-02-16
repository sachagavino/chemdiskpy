import glob, os, sys, shutil 
import numpy as np

from .. import radmc3d
from .. import nautilus
from .Grid import Grid

from .. constants.constants import autocm, M_sun, R_sun


class Model:

    def __init__(self):
        self.grid = Grid()

    def run_thermal(self, nphot=1e6, **keywords):
            self.run_thermal_radmc3d(nphot=nphot, **keywords)

    def run_localfield(self, nphot_mono=1e6):
        self.run_localfield_radmc3d(nphot_mono=nphot_mono)

    def run_thermal_radmc3d(self, nphot=1e6, verbose=True, timelimit=7200, \
            nice=None, **keywords):

        thermpath='thermal/'

        self.write_radmc3d(nphot_therm=nphot, **keywords)
        # ------- >  radmc3d.run.thermal(verbose=verbose, timelimit=timelimit, nice=nice)

        self.grid.temperature = radmc3d.read.dust_temperature(thermpath) # will be updated when multiple structures

        #-------- IF MULTIPLE STRUCTURES:
        # for i in range(len(self.grid.temperature)):
        #     nbspecies, n1, n2, n3 = self.grid.dustdensity[i].shape
        #     self.grid.temperature[i] = self.grid.temperature[i].reshape((nbspecies, n1, n2, n3))

        nbspecies, nr, nt, nph = self.grid.dustdensity[0].shape
        #self.grid.temperature[0] = self.grid.temperature[0].reshape((nbspecies, n1, n2, n3))
        self.grid.temperature[0] = np.reshape(self.grid.temperature[0], (nbspecies, nph, nt, nr))

    def run_localfield_radmc3d(self, nphot_mono=1e6, verbose=True, timelimit=7200):
        radmc3d.run.localfield(nphot_mono=nphot_mono, verbose=verbose, timelimit=timelimit)

        thermpath='thermal/'
        self.grid.localfield = radmc3d.read.localfield(thermpath)

        # #for i in range(len(self.grid.localfield)):
        # nbspecies, n1, n2, n3 = self.grid.dustdensity[0].shape
        # nlam = int(len(self.grid.localfield)/(n1*n2))
        # self.grid.localfield = np.reshape(radmc3d.read.localfield(thermpath), (nlam, n2, n1))
        # #self.grid.localfield = np.reshape(self.grid.localfield, (nlam, n2, n1, n3))
        # #print(self.grid.localfield.shape)
        # #self.grid.localfield.reshape((nlam, n1, n2, n3))


    def write_radmc3d(self, **keywords):
        print('writing RADMC3D input files...\n')
        #os.system("rm thermal/*.inp")
        if not os.path.exists('thermal'):
            os.makedirs('thermal')

        radmc3d.write.control(**keywords)

        mstar = []
        rstar = []
        xstar = []
        ystar = []
        zstar = []
        tstar = []

        
        for i in range(len(self.grid.stars)):
            mstar.append(self.grid.stars[i].mass*M_sun)
            rstar.append(self.grid.stars[i].radius*R_sun)
            xstar.append(self.grid.stars[i].x*autocm)
            ystar.append(self.grid.stars[i].y*autocm)
            zstar.append(self.grid.stars[i].z*autocm)
            tstar.append(self.grid.stars[i].temperature)

        radmc3d.write.stars(rstar, mstar, self.grid.lam, xstar, ystar, zstar, \
                tstar=tstar)

        radmc3d.write.wavelength_micron(self.grid.lam)

        radmc3d.write.mcmono_wavelength_micron(self.grid.monolam)

        if len(self.grid.isrf) != 0:
            radmc3d.write.external_rad(self.grid.isrf[0])
        
        if self.grid.coordsystem == 'spherical':
            radmc3d.write.amr_grid(self.grid.w1*autocm, self.grid.w2, self.grid.w3, gridstyle="regular", coordsystem=self.grid.coordsystem)

        radmc3d.write.dust_density(self.grid.dustdensity)

        dustopac = []
        filelist = glob.glob('thermal/dustkap*')
        for files in sorted(filelist):
            dustopac.append(files)
        radmc3d.write.dustopac(dustopac)

    # WRITE NAUTILUS INPUT FILES
    def write_nautilus(self, sizes=[0.1], uv_ref=3400, nH_to_AV_conversion=1.600e+21, rgrain=0.1, dtogas=1e-2, ref_radius=100,\
                       stop_time=3e6, static=True, param=True, element=True, abundances=True, \
                       activ_energies=True, surfaces=True,  network=True, grain_sizes=True,\
                       coupling_temp=True, coupling_av=True, **keywords):

        chemfold = 'chemistry'
        if os.path.exists(chemfold):
            shutil.rmtree(chemfold)
        os.makedirs(chemfold)

        nbspecies = len(sizes[-1])

        # coupling from RADMC3D output model toward NAUTILUS input model                
        if coupling_temp == True:
            T_dust = nautilus.coupling.dust_temperature(self.grid.temperature[0], self.grid.rchem*autocm, self.grid.zchem, self.grid.r*autocm, self.grid.theta, self.grid.hg_chem[0]) # dim(a, rchem, zchem)
        else:
            T_dust = np.expand_dims(self.grid.tgas_chem[0], axis=0) # expend to one extra dimension in order to match the shape of coupled T_dust.
        if coupling_av == True:
            av_z = nautilus.coupling.av_z(self.grid.localfield, self.grid.monolam, self.grid.rchem*autocm, self.grid.zchem, self.grid.r*autocm, self.grid.theta, self.grid.hg_chem[0]) # dim(rchem, zchem)
        else:
            av_z = self.grid.avz[0]

        print('writing NAUTILUS input files...\n')
        # write input NAUTILUS files
        for idx, r in enumerate(self.grid.rchem):
            path = 'chemistry' + '/' + str(r) + 'AU/'

            try:
                os.makedirs(path, exist_ok=True) 
            except IOError:
                print('There is no folder called chemistry.')
                sys.exit(1)
            uvflux = nautilus.write.uv_factor(uv_ref, ref_radius, r, self.grid.hg_chem[0][idx]/autocm)
            avnh_fact = nautilus.write.avnh_factor(nH_to_AV_conversion, dtogas, rgrain, self.grid.zchem)


            if network == True:
                nautilus.write.network(path)
            if element == True:
                nautilus.write.elements(path)
            if abundances == True:
                nautilus.write.abundances(path)
            if activ_energies == True:
                nautilus.write.activ_energies(path)
            if surfaces == True:
                nautilus.write.surfaces(path)

            if nbspecies == 1:
                if param == True:
                    nautilus.write.parameters(path, resolution=self.grid.nz_chem, stop_time=stop_time, uv_flux=uvflux, **keywords)
                if static == True:
                    nautilus.write.static(path, \
                                    self.grid.zchem, \
                                    self.grid.hg_chem[0][idx]/autocm, \
                                    self.grid.gasdensity_chem[0][idx,:], \
                                    self.grid.tgas_chem[0][idx], \
                                    av_z[idx, :], \
                                    T_dust[:,idx,:], \
                                    self.grid.dustdensity_chem[0][0,idx,:], \
                                    rgrain, \
                                    avnh_fact)
            if nbspecies > 1:
                if param == True:
                    nautilus.write.parameters_multi(path, resolution=self.grid.nz_chem, nb_grains=nbspecies, stop_time=stop_time, uv_flux=uvflux, **keywords)
                if static == True:
                    nautilus.write.static(path, \
                                    self.grid.zchem, \
                                    self.grid.hg_chem[0][idx]/autocm, \
                                    self.grid.gasdensity_chem[0][idx,:], \
                                    self.grid.tgas_chem[0][idx], \
                                    av_z[idx, :], \
                                    T_dust[:,idx,:], \
                                    self.grid.dustdensity_chem[0][0,idx,:], \
                                    rgrain, \
                                    avnh_fact)
                if grain_sizes == True:
                    nautilus.write.grain_sizes(path, sizes, self.grid.gasdensity_chem[0][idx,:], self.grid.dustdensity_chem[0][:,idx,:], T_dust[:,idx,:])


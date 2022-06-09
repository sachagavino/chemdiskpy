import glob, os, sys, shutil 
import numpy as np

from .. import radmc3d
from .. import nautilus
from .Grid import Grid

from .. constants.constants import autocm, M_sun, R_sun

import matplotlib.pyplot as plt



class Model:

    def __init__(self):
        self.grid = Grid()

    def thermal(self, nphot=1e6, run=True, write=True, **keywords):
        """ 
        Notes:
        run MC dust radiative transfer, open the resulting dust temperature as an array and computes the surface-area weigthed temperature. If run == False, user assumes the RADMC3D output files already exist.
        -----
	    """	
        if write == True:
            self.write_radmc3d(nphot_therm=nphot, **keywords)
        if run == True:
            self.run_thermal_radmc3d(nphot=nphot, **keywords)

            ##-------- IF MULTIPLE STRUCTURES:
            # for i in range(len(self.grid.temperature)):
            #     nbspecies, n1, n2, n3 = self.grid.dustdensity[i].shape
            #     self.grid.temperature[i] = self.grid.temperature[i].reshape((nbspecies, n1, n2, n3))

            
            ##WARNING: will not work for several structures. 
            # nb, nr, nt, nph = self.grid.dustdensity[0].shape
            # print(self.grid.dustdensity[0].shape)
            # nbspecies = 0
            # for istruct in range(len(self.grid.dustdensity)): 
            #     nbspecies += len(self.grid.dustdensity[istruct])

        if write == False:
            print('RADMC3D input files already present. Will continue.\n')

        thermpath='thermal/'    
        
        nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath)
        dust_density = radmc3d.read.dust_density(thermpath) # Gives a list of one numpy array. Will be updated when multiple structures
        nbspecies = int(len(dust_density[0])/(nx*ny*nz))

        dust_density[0] = np.reshape(dust_density[0], (nbspecies, nz, ny, nx))
        dust_density[0][dust_density[0]<=1e-100] = 1e-100
        self.grid.temperature = radmc3d.read.dust_temperature(thermpath) # Gives a list of one numpy array. Will be updated when multiple structures
        self.grid.temperature[0] = np.reshape(self.grid.temperature[0], (nbspecies, nz, ny, nx)) 

        dustmodel = self.grid.dust[0] # represents the first structure. If multiple structures, will not match with self.grid.temperature[0]
        dustmass = dustmodel.grainmass() # gram
        sizes = dustmodel.sizes()
        a = sizes[-1]*1e-4 # cm

        Ta_num = np.zeros((nz, ny, nx))
        Ta_denum = np.zeros((nz, ny, nx))

        for idx in range(len(a)):
            Ta_num += self.grid.temperature[0][idx, :, :, :]*(dust_density[0][idx, :, :, :]/dustmass[idx])*a[idx]**2
            Ta_denum += (dust_density[0][idx, :, :, :]/dustmass[idx])*a[idx]**2
        self.Ta = Ta_num/Ta_denum  #Ta is the area-weighted dust temperature.


    def localfield(self, nphot_mono=1e6, run=True):
        if run == True:
            self.run_localfield_radmc3d(nphot_mono=nphot_mono)

        thermpath='thermal/'
        self.grid.localfield = radmc3d.read.localfield(thermpath)
        nlam = len(self.grid.monolam)
        
        # #for i in range(len(self.grid.localfield)):
        nbspecies, nr, nt, nph = self.grid.dustdensity[0].shape
        # nlam = int(len(self.grid.localfield)/(nr*nt))
        #self.grid.localfield = np.reshape(radmc3d.read.localfield(thermpath), (nlam, n2, n1))
        self.grid.localfield = np.reshape(self.grid.localfield, (nlam, nph, nt, nr))
        # #self.grid.localfield.reshape((nlam, n1, n2, n3))

    def run_thermal_radmc3d(self, nphot=1e6, verbose=True, timelimit=7200, \
            nice=None, **keywords):
        radmc3d.run.thermal(verbose=verbose, timelimit=timelimit, nice=nice)


    def run_localfield_radmc3d(self, nphot_mono=1e6, verbose=True, timelimit=7200):
        radmc3d.run.localfield(nphot_mono=nphot_mono, verbose=verbose, timelimit=timelimit)


    def write_radmc3d(self, **keywords):
        print('\nWRITING RADMC3D INPUT FILES:')
        print('----------------------------\n')
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
    def write_nautilus(self, sizes=np.array([[0.1]]), uv_ref=3400, nH_to_AV_conversion=1.600e+21, rsingle=0.1, dtogas=1e-2, ref_radius=100,\
                       stop_time=3e6, nb_outputs = 64, static=True, param=True, element=True, abundances=True, \
                       activ_energies=True, surfaces=True, network=True, nmgc=False, single=False, \
                       coupling_temp=True, coupling_av=True, **keywords):

        chemfold = 'chemistry'
        if os.path.exists(chemfold):
            shutil.rmtree(chemfold)
        os.makedirs(chemfold)

        nbspecies = len(sizes[-1]) # if nbspecies = 1, then sizes = np.array([[self.rsingle]])
        # coupling from RADMC3D output model toward NAUTILUS input model                
        if coupling_temp == True:
            if not self.grid.temperature:
                print('The file thermal/dust_temperature.dat is not present or is corrupted.')
                sys.exit(1)
            elif self.grid.rchem.size == 0:
                print('The radii for the chemistry model (rchem) are not set. Please write something like m.grid.set_nautilus_grid(rchem, ...)')
                sys.exit(1)
            elif not self.grid.hg_chem:
                print('You did not define a NAUTILUS model. Please, write something like m.write_nautilus(sizes=sizes, ...).')
                sys.exit(1)
            T_dust = nautilus.coupling.dust_temperature(self.grid.temperature[0], self.grid.rchem*autocm, self.grid.zchem, self.grid.r*autocm, self.grid.theta, self.grid.hg_chem[0]) # dim(a, rchem, zchem)
            T_dust_single = nautilus.coupling.dust_temperature_single(self.Ta, self.grid.rchem*autocm, self.grid.zchem, self.grid.r*autocm, self.grid.theta, self.grid.hg_chem[0])
        else:
            T_dust = np.expand_dims(self.grid.tgas_chem[0], axis=0) # expend to one extra dimension in order to match the shape of coupled T_dust.

        if coupling_av == True:
            if not self.grid.localfield.size:
                print('The file thermal/mean_intensity.out is not present or is corrupted.')
                sys.exit(1)
            av_z = nautilus.coupling.av_z(self.grid.localfield, self.grid.monolam, self.grid.rchem*autocm, self.grid.zchem, self.grid.r*autocm, self.grid.theta, self.grid.hg_chem[0]) # dim(rchem, zchem)
        else:
            av_z = self.grid.avz[0]

        print('WRITING NAUTILUS INPUT FILES:')
        print('-----------------------------\n')
        # write input NAUTILUS files
        if nmgc == True:
             print('NMGC: Yes\n')
        else:
             print('NMGC: No\n')

        if single == True:
             print('Dust temperatures: area-weigthed\n')
        else:
             print('Dust temperatures: independent\n')

        for idx, r in enumerate(self.grid.rchem):
            path = 'chemistry' + '/' + str(r) + 'AU/'

            try:
                os.makedirs(path, exist_ok=True) 
            except IOError:
                print('There is no folder called chemistry.')
                sys.exit(1)

            uvflux = nautilus.write.uv_factor(uv_ref, ref_radius, r, self.grid.hg_chem[0][idx]/autocm)
            avnh_fact = nautilus.write.avnh_factor(nH_to_AV_conversion, dtogas, rsingle, self.grid.zchem)

            if nmgc == False:
                if param == True:
                    nautilus.write.parameters(path, nb_outputs=nb_outputs, resolution=self.grid.nz_chem, stop_time=stop_time, uv_flux=uvflux, **keywords)
                if static == True:
                    nautilus.write.static(path, \
                                    self.grid.zchem, \
                                    self.grid.hg_chem[0][idx]/autocm, \
                                    self.grid.gasdensity_chem[0][idx,:], \
                                    self.grid.tgas_chem[0][idx], \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], \
                                    self.grid.dustdensity_single_chem[0][idx,:], \
                                    rsingle, \
                                    avnh_fact)
            if nmgc == True:
                if param == True:
                    if nbspecies == 1 or single == True:
                        nautilus.write.parameters_nmgc(path, grain_temp='table_1D', nb_outputs=nb_outputs, resolution=self.grid.nz_chem, nb_grains=1, stop_time=stop_time, uv_flux=uvflux, **keywords)
                    else:
                        nautilus.write.parameters_nmgc(path, grain_temp='fixed_to_dust_size', nb_outputs=nb_outputs, resolution=self.grid.nz_chem, nb_grains=nbspecies, stop_time=stop_time, uv_flux=uvflux, **keywords)
                if static == True:
                    nautilus.write.static(path, \
                                    self.grid.zchem, \
                                    self.grid.hg_chem[0][idx]/autocm, \
                                    self.grid.gasdensity_chem[0][idx,:], \
                                    self.grid.tgas_chem[0][idx], \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], #if nbspecies > 1, the column is not read.\
                                    self.grid.dustdensity_single_chem[0][idx,:], #if nbspecies > 1, the column is not read. \
                                    rsingle, #if nbspecies > 1, the column is not read. \
                                    avnh_fact)
                if nbspecies > 1 and single == False:
                    #If single==True then 1D_grain_sizes.in is not read, but we create it any way. 
                    nautilus.write.grain_sizes(path, sizes, self.grid.gasdensity_chem[0][idx,:], self.grid.dustdensity_chem[0][:,idx,:], T_dust[:,idx,:])

            # if network == True:
            #     nautilus.write.network(path)
            if element == True:
                nautilus.write.elements(path)
            if abundances == True:
                nautilus.write.abundances(path)
            if activ_energies == True:
                nautilus.write.activ_energies(path)
            if surfaces == True:
                nautilus.write.surfaces(path)


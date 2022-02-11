import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust

m = modeling.YSOModel()  


#-----DISK PARAMETERS
nr, ntheta, nphi = 301, 181, 2
rin, rout= 1, 300
dtogas = 1e-2

#-----CREATE WAVELENGTH TABLE----
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) #microns
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) #microns


#-----CREATE RADMC GRID-----
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)

# #-----CREATE STARS-----
m.add_star(mass=0.5, luminosity=0.75, temperature=4000., x=0., y=0., z=0.)

# # #-----CREATE ISRF-----
#m.add_isrf(cut=2.e-1, d78=True, vdb82=True)

# # #-----CREATE DUST DENSITIES-----
d = dust.DustDistrib(rsingle=0.1, nb_sizes=1, rho_m=2.5) #microns
m.add_disk(dust=d, rin=rin, rout=rout, dtogas=dtogas, settling=True, coordsystem='spherical')
m.run_thermal(nphot = 1e6, \
              nphot_scat = 3e5, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 1, \
              tgas_eq_tdust = 0)

m.run_localfield(nphot_mono = 1000000)


rchem = [80, 100]
m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)

m.add_nautilusdisk(dust=d, dtogas=dtogas, settling=True)

m.write_nautilus(surface=False, abundances=False, parameters=False, uv_ref=3400, cr_ionisation_rate=1.3e-17, structure_type='1D_no_diff', dtogas=dtogas, rgrain=d.rsingle, ref_radius=m.disk.ref_radius, )

import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot

m = modeling.YSOModel()  


#-----DISK PARAMETERS
nr, ntheta, nphi = 301, 181, 2
rin, rout= 1, 500
dtogas = 1e-2
star_mass = 1.
dust_mass = 1e-4
q_c = 12

#-----CREATE WAVELENGTH TABLE----
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) #microns
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) #microns


#-----CREATE RADMC GRID-----
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)

# #-----CREATE STARS-----
m.add_star(mass=star_mass, luminosity=0.8, temperature=3900., x=0., y=0., z=0.)

# # #-----CREATE ISRF-----
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)

# -----CREATE DUST DENSITIES-----
d = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=2, rho_m=3.) #microns
sizes = d.sizes()
mass = d.grainmass()
print(sizes[-1])
m.add_disk(dust=d, rin=rin, rout=rout, star_mass=star_mass, dtogas=dtogas, dust_mass=dust_mass, q_c=q_c, settling=True, coordsystem='spherical')

# -----RUN THERMAL-----
m.thermal(run=True, nphot = 2e7, \
              nphot_scat = 2e7, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 2, \
              tgas_eq_tdust = 0)


m.localfield(run=True, nphot_mono = 2000000)

rchem = [10, 20, 30, 50, 80, 100, 140, 160, 200]
m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)

m.add_nautilusdisk(dust=d, rin=rin, rout=rout,  star_mass=star_mass, dust_mass=dust_mass, q_c=q_c, dtogas=dtogas, settling=True)

m.write_nautilus(sizes=sizes, uv_ref=3400, dtogas=dtogas, ref_radius=m.disk.ref_radius, cr_ionisation_rate=1.300E-17, network=False, coupling_av=False)

plot.vertical_temp(r=100)
plot.midplane_temp()
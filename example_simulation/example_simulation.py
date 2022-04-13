import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot

m = modeling.Structure()  


#-----DISK PARAMETERS
nr, ntheta, nphi = 301, 181, 2
rin, rout= 2, 500
dtogas = 1e-2
star_mass = 1.
dust_mass = 7.5e-5
q_c = 12
q_exp = 0.4

#-----CREATE WAVELENGTH TABLE----
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) #microns
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) #microns


#-----CREATE RADMC GRID-----
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)

# #-----CREATE STARS-----
m.add_star(mass=star_mass, luminosity=1.5, temperature=4100., x=0., y=0., z=0.)

# # #-----CREATE ISRF-----
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)

# # -----CREATE DUST DENSITIES-----
d = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=2, rho_m=3.) #microns
frac = d.massfraction()
sizes = d.sizes()
mass = d.grainmass()
print(sizes[-1])
print(frac)
print(mass)
m.add_disk(dust=d, rin=rin, rout=rout, Tmidplan_ref=15, star_mass=star_mass, q_exp= q_exp, dtogas=dtogas, dust_mass=dust_mass, q_c=q_c, settling=True, coordsystem='spherical')

# -----RUN THERMAL-----
m.thermal(run=True, nphot = 1e7, \
              nphot_scat = 1e7, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 2, \
              tgas_eq_tdust = 0)


m.localfield(run=False, nphot_mono = 2000000)

rchem = [50,55,60,100,110,120,130,140,150,160]
m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)

m.add_nautilusdisk(dust=d, rin=rin, rout=rout,  Tmidplan_ref = 15, star_mass=star_mass, dust_mass=dust_mass, q_c=q_c, q_exp= q_exp, dtogas=dtogas, settling=True)

m.write_nautilus(nmgc=True, sizes=sizes, uv_ref=3400, dtogas=dtogas, \
                 ref_radius=m.disk.ref_radius, is_h2_adhoc_form = 0, is_h2_formation_rate=1, \
                 height_h2formation=64, is_photodesorb = 1, cr_ionisation_rate=1.300E-17, \
                 stop_time=2e6, nb_outputs = 60, network=False, coupling_av=False)

plot.vertical_temp(r=100)
plot.density2D(d.grainmass())
plot.midplane_temp()
plot.temperature2D()
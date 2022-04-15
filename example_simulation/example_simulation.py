import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot

m = modeling.Structure()  


#-----PARAMETERS
nr, ntheta, nphi = 300, 301, 2
rin, rout= 1, 500
dtogas = 1e-2
star_mass = 1.5
disk_dust_mass = 1e-4 #total disk dust mass
env_dust_mass = 1e-5 #total envelope dust mass
q_c = 12
q_exp = 0.4

#-----CREATE WAVELENGTH TABLE----
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) #microns
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) #microns


#-----CREATE RADMC GRID-----
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)

# #-----CREATE STARS-----
m.add_star(mass=star_mass, luminosity=2, temperature=4000., x=0., y=0., z=0.)

# # #-----CREATE ISRF-----
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)

# # -----CREATE DUST DENSITIES-----
dust_disk = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=1, rho_m=3.) #microns
dust_envelope = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=1, rho_m=3.) #microns both structures has the same dust species.
frac = dust_disk.massfraction()
sizes = dust_disk.sizes()
mass = dust_disk.grainmass()


m.add_disk(dust=dust_disk, rin=rin, rout=rout, Tmidplan_ref=15, star_mass=star_mass, q_exp= q_exp, dtogas=dtogas, dust_mass=disk_dust_mass, q_c=q_c, settling=True, coordsystem='spherical')
m.add_envelope(dust=dust_envelope, rmin=1, rmax=1400, dust_mass=env_dust_mass, dtogas=dtogas, coordsystem='spherical')

# -----RUN THERMAL-----
m.thermal(run=True, nphot = 1e6, \
              nphot_scat = 1e6, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 2, \
              tgas_eq_tdust = 0)


# # #m.localfield(run=False, nphot_mono = 2000000)

# rchem = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360]
# m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)

# m.add_nautilusdisk(dust=d, rin=rin, rout=rout,  Tmidplan_ref = 15, star_mass=star_mass, dust_mass=dust_mass, q_c=q_c, q_exp= q_exp, dtogas=dtogas, settling=True)

# m.write_nautilus(nmgc=True, sizes=sizes, uv_ref=3400, dtogas=dtogas, \
#                  ref_radius=m.disk.ref_radius, is_h2_adhoc_form = 0, is_h2_formation_rate=1, \
#                  height_h2formation=64, is_photodesorb = 1, cr_ionisation_rate=1.300E-17, \
#                  stop_time=2e6, nb_outputs = 60, network=False, coupling_av=False)

# plot.vertical_temp(r=100)
#plot.density2D(d.grainmass())
plot.density2D(dust_disk.grainmass(), dust_envelope.grainmass())
#plot.midplane_temp()
plot.temperature2D()
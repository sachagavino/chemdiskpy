# chemdiskpy: Couple chemistry modeling with thermal modeling of protoplanetary disks with multiple grain sizes


This code is dedicated to help users create dust continuum radiative transfer simulations with RADMC3D and couple the output results to create input parameters for the gas-grain code NAUTILUS for protoplanetary disk models.
Although it is totally possible to use only one grain size and species, the main purpose of the package is to use multiple grain sizes, both for solving thermal and chemical problems.


## Installation

- last release from test [https://test.pypi.org/](https://test.pypi.org/):

        pip install -i https://test.pypi.org/simple/ chemdiskpy==0.1.0


## Quick start

### Set up a directory (this step is temporary and will not be required in the next version)
1. Create a working folder where you want to create a model and go to this folder.
2. TEMPORARY: Create manually a folder named **chemistry/** and another **thermal/**. This step is temporary and will be automatic in next version.
- **chemistry/** is where the NAUTILUS model is stored.
- **thermal/** is where the thermal model is stored.
3. In thermal folder, add an opacity table in the format of radmc3d using the script provided. and go back to the working directory.

## Set up a model
- All parameters can be written in the file parameters.py in the working directory. See folder example_simulation.

### IMPORT PACKAGE
- Open any notebook or script in the working directory.

```
import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot
```

### CREATE AN OBJECT:
```
m = modeling.YSOModel() 
```

### DISK PARAMETERS
```
nr, ntheta, nphi = 301, 181, 2 # number of points for the dust continuum radiative transfer
rin, rout= 1, 500 # in au
dtogas = 1e-2
star_mass = 0.8 # in solar mass 
```

### WAVELENGTH GRID
```
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) # in microns.
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) # in microns. This grid will be used for the computing of the local radiation field.
```

### RADMC3D GRID
- Typically, the grid is spherical for model with a central object like disks (but can be cartesian if needed):

```
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)
```

### STAR
- This is a simple black-body spectrum. UV from accreting material will be added in future versions. The star is located at coordinates (0,0,0)
```
m.add_star(mass=star_mass, luminosity=1., temperature=4500., x=0., y=0., z=0.)
```

### INTERSTELLAR RADIATION FIELD (optional)
- Draine 1978 between 91.2 and 200 nm, and with the extension of van Dishoeck & Black 1982 at longer wavelengths. Change 'cut' if you want another value.
```
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)
```

### DUST DISK STRUCTURE
- Start by creating a dust population. This library is originally made for multi-grain chemical code but this works fine with a single grain population.
- Set **nb_sizes=1** if you want to use one size. In that case, **amin** and **amax** will not be read and **rsingle** is the grain size.
```
d = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=2, rho_m=3.) #microns
sizes = d.sizes() # get the size 
mass = d.grainmass() # get the grain mass for each grain size. 
```

- Then, create a dust disk structure with spherical coordinates for RADMC3D with the values previously set. 
```
m.add_disk(dust=d, rin=rin, rout=rout, dtogas=dtogas, dust_mass=8e-5, settling=True, coordsystem='spherical')
```

### RUN THERMAL
- The main parameters are the number of photons **nphot** and **nphot_scat**. 
```
m.run_thermal(nphot = 5e5, \
              nphot_scat = 5e5, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 1, \
              tgas_eq_tdust = 0)
```

### PLOT THE RESULTS (optional)
- Dust temperature:
```
plot.temperature2D() 
```

- Dust density:
```
plot.density2D(mass) # Use the mass previously set up. 
```

- Radial dust temperature profile (in the midplane):
```
plot.midplane_temp() 
```

- Dust opacity (absorption, scattering, angles)
```
plot.opacity() 
```

### RUN LOCAL RADIATION (optional)
```
nphot_mono = 1e7 # choose number of photons
m.run_localfield(nphot_mono = nphot_mono)
```

### NAUTILUS GRID
- radii in au where you want to compute chemistry:
```
rchem = [10, 20, 30 , 40, 50, 60, 80, 100] 
```
- create NAUTILUS grid. 
- The value **max_h** is the maximum scale height (z/H) where chemistry will be computed. The grid is empty above. 
- The value **nz_chem** is the number of vertical points where chemistry is computed for each radius.
```  
m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)
```

### CREATE NAUTILUS DISK STRUCTURE
- Use the same dust model as for RADMC3D.
```
m.add_nautilusdisk(dust=d, dtogas=dtogas, settling=True)
```

### CREATE NAUTILUS DISK MODEL
- Create a ready-to-use nautilus disk model in **chemistry/**. The dust temperature (and local flux if needed) is extraced from the RADMC3D results.
```
m.write_nautilus(sizes=sizes, uv_ref=3400, dtogas=dtogas, ref_radius=m.disk.ref_radius, cr_ionisation_rate=1.300E-17, network=False, coupling_av=False)
```

- **network=False** implies that the user has to add his/her own chemical network. If True, the kida.14 network will be added in the user's model.
- WARNING: the parameter **coupling_av=True** works only if the local field is computed prior to this.

### PLOT THE RESULTS (optional)
- Vertical dust temperature:
```
plot.vertical_temp(r=100) 
```

### RUN NAUTILUS DISK MODEL
- Next step is to run as you would usually do the nautilus model in **chemistry/**.
- This will be added in a future update.
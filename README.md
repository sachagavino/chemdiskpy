# chemdiskpy: Couple chemistry modeling with thermal modeling of protoplanetary disks with multiple grain sizes


This code is dedicated to help users create dust continuum radiative transfer simulations with RADMC3D and couple the output results to create input parameters for the gas-grain code NAUTILUS for protoplanetary disk models.
Although it is totally possible to use only one grain size and species, the main purpose of the package is to use multiple grain sizes, both for solving thermal and chemical problems.


## Installation

- last release from test [https://test.pypi.org/](https://test.pypi.org/):

        pip install -i https://test.pypi.org/simple/ chemdiskpy==0.1.0


## Quick start

### Set up a directory (this step is temporary and will not be required in the next version)
- Create a working folder where you want to create a model and go to this folder.
- TEMPORARY: Create a folder named 'chemistry' and another 'thermal'. This step is temporary.
- in thermal folder, add a opacity table in the format of radmc3d using the script. and go back to the working directory.

## Set up a model

### IMPORT PACKAGE
Open any notebook or script in the working directory.

```
import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot
```

### Create a object:
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

### RADMC GRID
Typically, the grid is spherical for model with a central object like disks:

```
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)
```

### STAR
This is a simple black-body spectrum. UV from accreting material will be added in future versions. The star is located at coordinates (0,0,0)
```
m.add_star(mass=star_mass, luminosity=1., temperature=4500., x=0., y=0., z=0.)
```

### INTERSTELLAR RADIATION FIELD
Draine 1978 between 91.2 and 200 nm, and with the extension of van Dishoeck & Black 1982 at longer wavelengths. Change 'cut' if you want another value.
```
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)
```

### DUST DISK STRUCTURE
```
Start with a dust populations. This library is originaly made for multi-grain chemistry code but this works fine with a single grain population.
Choose ```nb_sizes=1``` if you want to use one size. In that case amin and amax will not be read. ```rsingle``` is the grain size.
d = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=2, rho_m=3.) #microns
sizes = d.sizes() # get the size 
mass = d.grainmass() # get the grain mass for each grain size. 

Then, create a dust disk structure with spherical coordinates for RADMC3D with the values previously set. 
m.add_disk(dust=d, rin=rin, rout=rout, dtogas=dtogas, dust_mass=8e-5, settling=True, coordsystem='spherical')
```


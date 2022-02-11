# chemdiskpy: Couple chemistry modeling with thermal modeling of protoplanetary disks with multiple grain sizes


This code is dedicated to help users create dust continuum radiative transfer simulations with RADMC3D and couple the output results to create input parameters for the gas-grain code NAUTILUS for protoplanetary disk models.
Although it is totally possible to use only one grain size and species, the main purpose of the package is to use multiple grain sizes, both for solving thermal and chemical problems.


## Installation

- last release from test [https://test.pypi.org/](https://test.pypi.org/):

        pip install -i https://test.pypi.org/simple/ chemdiskpy==0.1.0


## Quick start
- Set up a model:

##### Import chemdiskpy packages
```
import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot
```

##### Create a object:
```
m = modeling.YSOModel() 
```

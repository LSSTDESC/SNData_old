# SNData
A repository with data structures useful for SN analysis and simulations.

## Description and Rationale
The basic input to most (if not all) analysis codes in supernovae are light curves and metadata for individual supernovae. These can be collected from real data or simulations, and each dataset/simulation might have its own format. The objectives here are
- to be able to interact with some inputs and obtain light curves and metadata as appropriate for our analysis codes
- to be able to write out simulations in a form that is useful for our simulation codes 

### Guiding Rationale
- We should plan for being able to use datasets of the size expected for LSST, without requiring everyone to have to do this from the start or when using this to analyze older datasets. For simulated datasets of the size of LSST, the answer is clearly a set of database tables or HDF5 files for portable, appropriate subsamples.
- We should plan to be able to interact with the LSST DM databases as a source of light curves. This is being developed as part of the [LSST DESC Monitor](https://github.com/LSSTDESC/Monitor) package, and we should be able to interface with this.

### Functionality

#### Read
In terms of simulation formats that we would like to be able to read:
- [SNANA](http://snana.uchicago.edu/) fits and ASCII formats
- SALT format (for example the JLA dataset)
- SNsims hdf5 files (and database formats when they become available)

#### API for light curves:
- Light Curves for basic data
- Light Curves for coadded data
- SNCosmo format light curves (ie. `astropy.table.Table`) available for fits 
#### Write
 We will also use these formats to write datasets for simulations.

## Requirements
- throughputs : Either install and setup the sims stack 
## Installation

## Contributing :
To contribute to the code:
1. please create a github issue describing the changes to make 

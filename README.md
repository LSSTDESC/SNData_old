# SNData
A repository with data structures useful for SN analysis and simulations.

## Installation
See the [install file](./install/installation.md)
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

In order to accept light curves from some data source, we would need the following minimal set of columns. Additional information is welcome, but these are the minimal number of columns:

_Light Curves_ :

In order to be able to connect to DM processed light curves, the light curves could be read in the format of an astropy Table with a minimal set of columns, and a metadata dictionary. 

The columns are described below, and also the keys to the dicttionary: 
band, mjd, flux, flux_err, zp, zpsys 

If calibrated zp can be et to zero. 
In case of LSST, a visit ID should be provided additionally (or in place of) band.

The astropy table should come with the following minimal metadata
snid, ra, dec, processing_tags 

 
#### Write
 We will also use these formats to write datasets for simulations.

## Requirements
- throughputs : Either install and setup the sims stack  
## Contributing :
To contribute to the code:
1. please create a github issue describing the changes to make 

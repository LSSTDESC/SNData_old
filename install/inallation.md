In order to install SNData

1. First clone the lsst throughputs package within a directory `mydir`. `mydir` is the absolute path to a location where the THROUGHPUTS are being cloned.
```
cd mydir
git clone git@github.com:lsst/throughputs.git
export THROUGHPUTS_DIR='mydir/THROUGHPUTS'
```

2. clone this directory and run `python setup.py install --user`


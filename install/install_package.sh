#!/usr/bin/env bash
pushd ../
git clone  https://github.com/lsst/throughputs
cd throughputs
export THROUGHPUTS_DIR=`pwd`
echo "THROUGHPUTS set to " $THROUGHPUTS_DIR
popd 
python setup/generate_requirements.py
python setup.py install --user

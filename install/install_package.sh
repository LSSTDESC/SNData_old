#!/usr/bin/env bash
pushd ../
git clone  https://github.com/lsst/throughputs
cd throughputs
export THROUGHPUTS_DIR=`pwd`
popd 
python setup/generate_requirements.py
python setup.py install --user

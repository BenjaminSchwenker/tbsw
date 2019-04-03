#!/bin/bash



if [ -z "$ROOTSYS" ]
then
      echo "\$ROOTSYS is not set (or empty). Please source ROOTSYS"
      exit -1
fi

# make sure we stop immediately on errors
set -e
 

mkdir -p build
pushd build
echo "Setting up build: >>cmake ..<<"
cmake ..
echo "Building: >>make -j8<<"
make -j8
popd

echo "Copying init script: workspace/init_tbsw.sh "
cp build/init_tbsw.sh workspace/init_tbsw.sh
source workspace/init_tbsw.sh

############################################################
# Create some workspace folders with special meaning in tbsw

# Folder to collect all calibration tags
mkdir -p workspace/localDB

# Folder to put logs for all processed runs
mkdir -p workspace/tmp-runs

# Folder to collect reconstructed DUT root files 
mkdir -p workspace/root-files

# Folder to collect DQM pdfs 
mkdir -p workspace/results


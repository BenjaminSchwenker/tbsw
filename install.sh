#!/bin/bash



if [ -z "$ROOTSYS" ]
then
      echo "\$ROOTSYS is not set (or empty). Please source Rootsys"
      return -1
fi

# make sure we stop immediately on errors
set -e
 

if [ -d "build" ]; then
  echo "Remove build dir"
  #rm -fr build/* 
fi
mkdir -p build
pushd build
echo "Setting up build: >>cmake ..<<"
cmake ..
echo "Building: >>make -j8<<"
make -j3
popd

echo "Creating init script: workspace/init_tbsw.sh "

############################################################
# Please do not edit stuff below here !!!!!!!!!!!!!!!!!!!!!!
export TBSW_HOME=$PWD/build
export TBSW_BIN=$TBSW_HOME/bin 
export TBSW_LIB=$TBSW_HOME/lib  
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH


echo "#!/bin/bash" > workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "#    ROOT                                                                        " >> workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "export ROOTSYS="${ROOTSYS}"" >> workspace/init_tbsw.sh
echo "export PATH="${ROOTSYS}/bin:\$PATH"" >> workspace/init_tbsw.sh
echo "export LD_LIBRARY_PATH="${ROOTSYS}/lib:\$LD_LIBRARY_PATH"" >> workspace/init_tbsw.sh
echo "export PYTHONPATH="${ROOTSYS}/lib:$PWD/source:\$PYTHONPATH"" >> workspace/init_tbsw.sh
echo "export ROOT_INCLUDE_PATH="${TBSW_LIB}"" >> workspace/init_tbsw.sh
echo "" >> workspace/init_tbsw.sh

echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "#    TBSW                                                                        " >> workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "export PATH="$TBSW_HOME/bin:\$PATH"" >> workspace/init_tbsw.sh
echo "export LD_LIBRARY_PATH="$TBSW_HOME/lib:\$LD_LIBRARY_PATH"" >> workspace/init_tbsw.sh
echo "export MARLIN_DLL="${TBSW_HOME}/lib/libTBReco.so:${TBSW_HOME}/lib/libEudaqInput.so:"" >> workspace/init_tbsw.sh
echo "export MARLIN="${TBSW_HOME}"" >> workspace/init_tbsw.sh 
echo "" >> workspace/init_tbsw.sh


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


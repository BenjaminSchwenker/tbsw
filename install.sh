#!/bin/bash

export CLHEP_HOME=/home/benjamin/work/CLHEP-2-3-4-5/lib/CLHEP-2.3.4.5
export ROOTSYS=/home/benjamin/work/root-v6-10-08

############################################################
# Please do not edit stuff below here !!!!!!!!!!!!!!!!!!!!!!
export TBSW_HOME=$PWD/source  
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

############################################################
# Install gear
cd source/gear
mkdir build
cd build
cmake ..
make install -j4
cd ../..

############################################################
# Install lcio
cd lcio
mkdir build
cd build 
cmake ..
make install  -j4
cd ../..

############################################################
# Install Marlin  
cd Marlin
mkdir build
cd build 
cmake -C ../BuildSetup.cmake ..
make install -j4
cd ../.. 

############################################################
# Install EudaqInput
cd EudaqInput
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install  -j4
cd ../..

############################################################
# Install TBTools  
cd TBTools
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install -j4
cd ../..

############################################################
# Install TBReco  
cd TBReco
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install -j4
cd ../../..         


############################################################
# Create bash script to set env variables  
echo "#!/bin/bash" > workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "#    ROOT                                                                        " >> workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "export ROOTSYS="${ROOTSYS}"" >> workspace/init_tbsw.sh
echo "export PATH="${ROOTSYS}/bin:${PATH}"" >> workspace/init_tbsw.sh
echo "export LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"" >> workspace/init_tbsw.sh
echo "export PYTHONPATH="${ROOTSYS}/lib"" >> workspace/init_tbsw.sh
echo "" >> workspace/init_tbsw.sh

echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "#    LCIO                                                                        " >> workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "export LCIO="${TBSW_HOME}/lcio"" >> workspace/init_tbsw.sh
echo "export PATH="$LCIO/tools:$LCIO/bin:$PATH"" >> workspace/init_tbsw.sh
echo "export LD_LIBRARY_PATH="$LCIO/lib:$LD_LIBRARY_PATH"" >> workspace/init_tbsw.sh
echo "" >> workspace/init_tbsw.sh

echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "#    Marlin                                                                      " >> workspace/init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> workspace/init_tbsw.sh
echo "export MARLIN="${TBSW_HOME}/Marlin"" >> workspace/init_tbsw.sh
echo "export PATH="${MARLIN}/bin:${PATH}"" >> workspace/init_tbsw.sh
echo "export MARLIN_DLL="${TBSW_HOME}/TBReco/lib/libTBReco.so:${TBSW_HOME}/EudaqInput/lib/libEudaqInput.so:"" >> workspace/init_tbsw.sh
echo "" >> workspace/init_tbsw.sh

############################################################
# Create some workspace folders with special meaning in tbsw

# Folder to collect all calibration tags
mkdir workspace/localDB

# Folder to put logs for all processed runs
mkdir workspace/tmp-runs

# Folder to collect reconstructed DUT root files 
mkdir workspace/root-files

# Folder to collect DQM pdfs 
mkdir workspace/results


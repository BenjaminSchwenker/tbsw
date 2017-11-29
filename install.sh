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
make install
cd ../..

############################################################
# Install lcio
cd lcio
mkdir build
cd build 
cmake ..
make install 
cd ../..

############################################################
# Install Marlin  
cd Marlin
mkdir build
cd build 
cmake -C ../BuildSetup.cmake ..
make install
cd ../.. 

############################################################
# Install EudaqInput
cd EudaqInput
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install 
cd ../..

############################################################
# Install TBTools  
cd TBTools
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install
cd ../..

############################################################
# Install TBReco  
cd TBReco
mkdir build
cd build 
cmake -C ../config/BuildSetup.cmake ..
make install
cd ../../..         

############################################################
# Create bash script to set env variables  

if [ -f "init_tbsw.sh" ]
then
	rm init_tbsw.sh
fi

echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "#    ROOT                                                                        " >> init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "export ROOTSYS="${ROOTSYS}"" >> init_tbsw.sh
echo "export PATH="${ROOTSYS}/bin:${PATH}"" >> init_tbsw.sh
echo "export LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"" >> init_tbsw.sh
echo "export PYTHONPATH="${ROOTSYS}/lib"" >> init_tbsw.sh
echo "" >> init_tbsw.sh

echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "#    LCIO                                                                        " >> init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "export LCIO="${TBSW_HOME}/lcio"" >> init_tbsw.sh
echo "export PATH="$LCIO/tools:$LCIO/bin:$PATH"" >> init_tbsw.sh
echo "export LD_LIBRARY_PATH="$LCIO/lib:$LD_LIBRARY_PATH"" >> init_tbsw.sh
echo "" >> init_tbsw.sh

echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "#    Marlin                                                                      " >> init_tbsw.sh
echo "#--------------------------------------------------------------------------------" >> init_tbsw.sh
echo "export MARLIN="${TBSW_HOME}/Marlin"" >> init_tbsw.sh
echo "export PATH="${MARLIN}/bin:${PATH}"" >> init_tbsw.sh
echo "export MARLIN_DLL="${TBSW_HOME}/TBReco/lib/libTBReco.so:${TBSW_HOME}/EudaqInput/lib/libEudaqInput.so:"" >> init_tbsw.sh
echo "" >> init_tbsw.sh


cp init_tbsw.sh workspace



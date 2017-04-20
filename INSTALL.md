
### INSTALL test beam analysis software (TBSW) ###

The installation procedure is tested on ubuntu 16.04. 


# Prerequisites #

Make sure the following packages are installed on your local machine before you try to install TBSW. 

# Cmake #

You need version of Cmake greater 3. You can find more information about Cmake at this URL https://cmake.org/download/. 

In case you are working with Ubuntu, the following recipe may work for you as well: 

- sudo apt-get install software-properties-common

- sudo add-apt-repository ppa:george-edison55/cmake-3.x

- sudo apt-get update

- sudo apt-get install cmake

taken from http://askubuntu.com/questions/610291/how-to-install-cmake-3-2-on-ubuntu-14-04

# CLHEP #

The project home page can be found at this URL http://proj-clhep.web.cern.ch/proj-clhep/. It is most convinient to install CLHEP from the git 
repository. 

- git clone https://gitlab.cern.ch/CLHEP/CLHEP.git

- mkdir CLHEP_Install

- mkdir CLHEP_Build

- cd CLHEP_Build

- cmake -DCMAKE_INSTALL_PREFIX=/home/fmu/CLHEP_Install /home/fmu/CLHEP

- cmake --build . --config RelWithDebInfo

- ctest

- cmake --build . --target install

# JAVA #

sudo apt-get install default-jdk

# Root #

The project home page can be found at this URL https://root.cern.ch/downloading-root. One approach to get the source code is to use the public GIT
repository to get the latest version.

- git clone http://root.cern.ch/git/root.git

The release specific tag can be obtained using for example:

- cd root
- git checkout -b v6-08-06 v6-08-06  

You want to install in a generic directory, depending on environment variables ROOTSYS, LD_LIBRARY_PATH, and PATH.

- mkdir <builddir>
- cd <builddir>
- cmake ../root
- cmake --build . 

Add bin/ to PATH and lib/ to LD_LIBRARY_PATH. For the sh shell family do:
   
- . bin/thisroot.sh

# Python #

There are many ways to get python. One convinient way is to use the ANACONDA package that can be found here https://www.continuum.io/DOWNLOADS. 

# QT4 #

sudo apt install libqt4-dev    

# TBSW #  
 
You can get the TBSW source code from a puplic git repository at bitbucket.

- git clone https://BenjaminSchwenker@bitbucket.org/BenjaminSchwenker/tbsw.git

- cd tbsw


Now, open the script install.sh and edit the two first lines with exports for ROOTSYS and CLHEP to the locations on your local machine. Save the edits and 
close the file. Run the install script: 

- . install.sh

After the installation is done, you can test the software by running a small test beam simulation. 

- cp -r workspace /home/myself/workspace-test

- cd /home/me/workspace-test 

- . init_tbsw.sh 

All actuall work with tbsw is encapsulated in workspaces. The installation comes with a workspace template workspace. The script init_tbsw.sh sets all needed
environment variables needed for working with TBSW.

- ./simulate.py -g steering-files/x0-sim/gear-misaligned.xml -x steering-files/x0-sim/simulation.xml -o mctest.slcio

- ./calibrate.py -i mctest.slcio -x steering-files/x0-sim/ -c test

- ./reco.py -i mctest.slcio -x steering-files/x0-sim/reco.xml -c test

The first command creates a file mctest.slcio containing a number of events with tracks crossing a pixel telescope. The second command calibrates the telescope 
from data and creates a calibration tag 'test'. The thrid command uses the calibration tag 'test' to reconstruct all simulated events. The result is a root file 
in the folder root-files.  

Some more information on the usage of tbsw can be found in the README files. 


Thats all;)	

benjamin.schwenker@phys.uni-goettingen.de











4) TBSW:

git clone https://fmu984@bitbucket.org/BenjaminSchwenker/tbsw.git

install.sh

edit: export CLHEP_HOME=/home/fmu/basf2/externals/v01-01-07/Linux_x86_64/opt/lib/CLHEP-2.2.0.4

edit: export ROOTSYS=/home/fmu/root_v5.34.24

cd tbsw cd workspace source init_tbsw.sh

./calibrate.py -i /home/fmu/TB_Benjamin/data/run000096_svdpxdteldigits_35k.slcio -x steering-files/vxd-april16-run96

./reco.py -i '/home/fmu/TB_Benjamin/data/run000096_svdpxdteldigits_35k.slcio' -x 'steering-files/vxd-april16-run96/reco.xml'


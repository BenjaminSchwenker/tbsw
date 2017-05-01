
# INSTALL test beam analysis software (TBSW) 

The installation procedure is tested on ubuntu 16.04. 


# Prerequisites 

Make sure the following packages are installed on your local machine before you try to install TBSW. 

# Cmake 

You need version of Cmake greater 3. You can find more information about Cmake at this URL https://cmake.org/download/. 

In case you are working with Ubuntu, the following recipe may work for you as well: 

```
$ sudo apt-get install software-properties-common
```
```
$ sudo add-apt-repository ppa:george-edison55/cmake-3.x
```
```
$ sudo apt-get update
```
```
$ sudo apt-get install cmake
```

taken from http://askubuntu.com/questions/610291/how-to-install-cmake-3-2-on-ubuntu-14-04

# CLHEP 

The project home page can be found at this URL http://proj-clhep.web.cern.ch/proj-clhep/. It is most convinient to install CLHEP from the git 
repository. 

```
$ git clone https://gitlab.cern.ch/CLHEP/CLHEP.git
$ mkdir CLHEP_Install
$ mkdir CLHEP_Build
$ cd CLHEP_Build
$ cmake -DCMAKE_INSTALL_PREFIX=/home/fmu/CLHEP_Install /home/fmu/CLHEP
$ cmake --build . --config RelWithDebInfo
$ ctest
$ cmake --build . --target install
```

# JAVA #

- sudo apt-get install default-jdk

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

- sudo apt install libqt4-dev    

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

All actuall work with tbsw is encapsulated in workspaces. The installation comes with a template workspace. You can test the installation by entering the 
workspace folder and do: 

- . init_tbsw.sh
- python tbsw_example.py


The script init_tbsw.sh sets all needed environment variables needed for working with TBSW. The script tbsw_example.py shows how to work with TBSW using 
python scripting. 


Some more information on the usage of tbsw can be found in the README files. 


Thats all;)	

benjamin.schwenker@phys.uni-goettingen.de



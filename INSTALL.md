
# INSTALL test beam software (tbsw) 

The installation procedure is tested on ubuntu 16.04. The requirements for an installation on Scientific Linux 7 are explained below.


# Prerequisites 

Make sure the following packages are installed on your local machine before you try to install tbsw. 

# Cmake 

You need version of Cmake greater 3. You can find more information about Cmake at this URL https://cmake.org/download/. 

In case you are working with Ubuntu, the following recipe may work for you as well: 

```
sudo apt install cmake
```

# Python 

Create a conda environment with a recent version of python3 e.g. 3.9:
```
conda create --name tbsw python=3.9
```
Then install the numpy and scipy packages for python:
```
conda install numpy scipy
```

# Root 

Install a recent version of root from looking at the project home page URL https://root.cern.ch/downloading-root. The release 6.16/00 was tested to be working fine. 

You want to install root in a generic directory, depending on environment variables ROOTSYS, LD_LIBRARY_PATH, and PATH.

1) Get the sources of the latest ROOT (see above). This folder will be called `./root` here. 

2) Type the build commands:
```
mkdir <builddir>
cd <builddir>
cmake ../root
cmake --build . [ -- -j<N> ] [ or simply "make -j<N>" on Unix systems ] 
```

3) Add bin/ to PATH and lib/ to LD_LIBRARY_PATH. For the sh shell family do:
```
. bin/thisroot.sh
```

This script also exports the ROOTSYS environment variable needed later for installing tbsw. 

4) try running root:
```
root
```

# TBSW   
 
You can get the tbsw source code from a puplic git repository at bitbucket.

```
git clone https://github.com/BenjaminSchwenker/tbsw.git
cd tbsw
```

Please make sure that the path to your ROOT instalation is set in $ROOTSYS
Run the install script:

```
./install.sh
```

This script will create a build folder "build" and will compile tbsw.
If you want to clean the tbsw directory for a complete rebuild, you can clean the tbsw directory with:

```
./make_clean.sh
```

Alternative:
If you want to build TBSW in a different folder, want a special build type or a special Root version to be used, please use:

```
cd your/build/dir
cmake path/to/tbsw -DROOT_HOME=/your/root/folder -DCMAKE_INSTALL_PREFIX=Release
cmake --build .
```

After the installation is done, you can test the software by running a small test beam simulation. All actuall work with tbsw is encapsulated in workspaces. The installation comes with a template workspace. The first two lines make a copy of the template workspace and cd into it. 

```
cp -r workspace ~/workspace-test
cd ~/workspace-test 
```

The script init_tbsw.sh sets all needed environment variables. The script tbsw_example.py shows how to run the simulation using python scripting. 

```
. init_tbsw.sh 
python example.py
```

# INSTALL test beam software on Scientific Linux 7

For Scientific Linux 7 only GCC 4.8 is installed, but tbsw uses c++14 features (lambda function parameter type capture with auto), so a newer version is needed. 
It can be obtained with the devtoolset:

```
# Install Software Collections in Scientific Linux 7
sudo yum install yum-conf-repos
sudo yum install yum-conf-softwarecollections

# Install Developer Toolset 7
sudo yum install devtoolset-7

# Enter Developer Toolset 7 Environment
scl enable devtoolset-7 bash

# Check your gcc version
gcc --version
```
Afterwards just follow the instructions for the installation on Ubuntu.


Some more information on the usage of tbsw can be found in the README files. 


Thats all;)	

benjamin.schwenker@phys.uni-goettingen.de



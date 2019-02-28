
# INSTALL test beam software (tbsw) 

The installation procedure is tested on ubuntu 16.04. The requirements for an installation on Scientific Linux 7 are explained below.


# Prerequisites 

Make sure the following packages are installed on your local machine before you try to install TBSW. 

# Cmake 

You need version of Cmake greater 3. You can find more information about Cmake at this URL https://cmake.org/download/. 

In case you are working with Ubuntu, the following recipe may work for you as well: 

```
$ sudo apt install cmake
```

# Python 

You should use your default system version of Python. Check that library python-dev is installed on your system to work with PyRoot. You can install by typing: 

```
sudo apt-get install python-dev
```

A few more python modules are needed to run the example scripts: 

```
sudo apt-get install python-numpy
sudo apt-get install python-scipy
```

# Root 

The project home page can be found at this URL https://root.cern.ch/downloading-root. One approach to get the source code is to use the public GIT
repository to get the latest version.

```
$ git clone http://root.cern.ch/git/root.git
```

The release specific tag can be obtained using for example:

```
$ cd root
$ git checkout -b v6-17-01 v6-17-01  
$ cd ..
```

You want to install in a generic directory, depending on environment variables ROOTSYS, LD_LIBRARY_PATH, and PATH.

```
$ mkdir <root-installdir>
$ mkdir <root-builddir>
$ cd <root-builddir>
$ cmake -DCMAKE_INSTALL_PREFIX=<root-installdir>  ../root
$ cmake --build . --target install
```

The folder <root-builddir> may be deleted afterwards. In the end, the folder <root-installdir> should contain bin/ include/ and lib/ subfolders.
In order to test the root installation, just do:

 
```  
$ cd <absolute-path-to-root-installdir>
$ . bin/thisroot.sh
$ root 
```

# TBSW   
 
You can get the TBSW source code from a puplic git repository at bitbucket.

```
$ git clone https://BenjaminSchwenker@bitbucket.org/BenjaminSchwenker/tbsw.git
$ cd tbsw
```

Please make sure that the path to your ROOT instalation is set in $ROOTSYS
Run the install script:

```
$ . install.sh
```

This script will create a build folder "build" and will compile tbsw.
If you want to clean the tbsw directory for a complete rebuild, you can clean the tbsw directory with:

```
$ . make_clean.sh
```

Alternative:
If you want to build TBSW in a different folder, want a special build type or a special Root version to be used, please use:

```
$ cd your/build/dir
$ cmake path/to/tbsw -DROOT_HOME=/your/root/folder -DCMAKE_INSTALL_PREFIX=Release
$ cmake --build .
```

After the installation is done, you can test the software by running a small test beam simulation. All actuall work with tbsw is encapsulated in workspaces. The installation comes with a template workspace. The first two lines make a copy of the template workspace and cd into it. 

```
$ cp -r workspace ~/workspace-test
$ cd ~/workspace-test 
```

The script init_tbsw.sh sets all needed environment variables. The script tbsw_example.py shows how to run the simulation using python scripting. 

```
$ . init_tbsw.sh 
$ python tbsw_example.py
```

# INSTALL test beam software on Scientific Linux 7

For Scientific Linux 7 only GCC 4.8 is installed, but tbsw uses c++14 features (lambda function parameter type capture with auto), so a newer version is needed. 
It can be obtained with the devtoolset:

```
# Install Software Collections in Scientific Linux 7
$ sudo yum install yum-conf-repos
$ sudo yum install yum-conf-softwarecollections

# Install Developer Toolset 7
$ sudo yum install devtoolset-7

# Enter Developer Toolset 7 Environment
$ scl enable devtoolset-7 bash

# Check your gcc version
$ gcc --version
```
Afterwards just follow the instructions for the installation on Ubuntu.


In case the script is working, your console should display the following: 

```
[INFO] Starting to simulate mc.slcio ...
[INFO] Marlin simulation.xml is done
[INFO] Starting to calibrate file mc.slcio ...
[INFO] Marlin hotpixelkiller.xml is done
[INFO] Marlin cluster-calibration-mc.xml is done
[INFO] Marlin correlator.xml is done
[INFO] Marlin kalmanalign-iteration-1.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin telescope-dqm.xml is done
[INFO] Created new caltag  default
[INFO] Done processing file mc.slcio
[INFO] Starting to reconstruct file mc.slcio ...
[INFO] Marlin reco.xml is done
[INFO] Done processing file mc.slcio
```

Some more information on the usage of tbsw can be found in the README files. 


Thats all;)	

benjamin.schwenker@phys.uni-goettingen.de



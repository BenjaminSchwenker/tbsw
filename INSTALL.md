
# INSTALL test beam software (tbsw) 

The installation procedure is tested on ubuntu 16.04. 


# Prerequisites 

Make sure the following packages are installed on your local machine before you try to install TBSW. 

# Cmake 

You need version of Cmake greater 3. You can find more information about Cmake at this URL https://cmake.org/download/. 

In case you are working with Ubuntu, the following recipe may work for you as well: 

```
$ sudo apt install cmake
```


# CLHEP 

The project home page can be found at this URL http://proj-clhep.web.cern.ch/proj-clhep/. It is most convinient to install CLHEP from the git 
repository. 

```
$ git clone https://gitlab.cern.ch/CLHEP/CLHEP.git
```

This create a folder CLHEP in your current working directory containing the CLHEP repository. Now create two further directories <installdir> and <builddir>
and build the software.  


```
$ mkdir <clhep-installdir> 
$ mkdir <clhep-builddir> 
$ cd <clhep-builddir> 
$ cmake -DCMAKE_INSTALL_PREFIX=<clhep-installdir>  ../CLHEP
$ cmake --build . 
$ ctest
$ cmake --build . --target install
```

The folder <clhep-builddir> may be deleted afterwards. In the end, the folder <clhep-installdir> should contain bin/ include/ and lib/ subfolders. 

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
$ git checkout -b v6-10-08 v6-10-08  
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

Open the script install.sh and edit the two first lines with exports for ROOTSYS and CLHEP_HOME to the locations on your local machine. 
For example: 

```
export CLHEP_HOME=<absolute-path-to-clhep-installdir>/lib/CLHEP-2.3.4.3         # See CLHEP section above, points to folder containing CLHEPConfig.cmake
export ROOTSYS=<absolute-path-to-root-installdir>                               # See Root section above
```

Save the edits and close the file. Run the install script: 

```
$ . install.sh
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



#############################################################################
# cmake example build setup for Marlin
# edit accordingly - pathes to packages !
# 
#
# For building Marlin with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
#         or to use an existing installation
#       cmake -C ../BuildSetup.cmake -C $ILCRELEASE/vxx-yy-zz/ILCSoft.cmake ..
# (4) $ make install
#
# @author Jan Engels, DESY
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

# Path to LCIO
SET( LCIO_HOME "$ENV{TBSW_HOME}/lcio"
    CACHE PATH "Path to LCIO" FORCE )

# Path to GEAR
SET( GEAR_HOME "$ENV{TBSW_HOME}/gear"
    CACHE PATH "Path to GEAR" FORCE )

# Path to LCCD
SET( LCCD_HOME "$ENV{TBSW_HOME}/lccd"
    CACHE PATH "Path to LCCD" FORCE )

# CMake Modules Path
SET( CMAKE_MODULE_PATH "$ENV{TBSW_HOME}/CMakeModules"
    CACHE PATH "Path to CMake Modules" FORCE )

#############################################################################
# Marlin GUI
#############################################################################

# You need to set QT4_HOME path if you want to build the GUI
SET( MARLIN_GUI ON CACHE BOOL "Set to ON to build Marlin GUI" FORCE )


#############################################################################
# Project options
#############################################################################

#SET( INSTALL_DOC OFF CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )

# set cmake build type, default value is: RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
#SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE )

#############################################################################
# Advanced options
#############################################################################

# You can only a static Marlin binary under Linux!!!
#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

#SET( MARLIN_NO_DLL ON CACHE BOOL "Set to ON to build Marlin without DLL support" FORCE )

# installation path for Marlin
#SET( CMAKE_INSTALL_PREFIX "/foo/bar" CACHE STRING "Where to install Marlin" FORCE )

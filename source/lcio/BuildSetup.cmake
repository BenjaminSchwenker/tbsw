#############################################################################
# cmake build setup for LCIO
#
# For building lcio with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
# @author Jan Engels, DESY
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

SET( BUILD_32BIT_COMPATIBLE OFF 
    CACHE BOOL "Set to Off for 64bit build" FORCE   )
 
#############################################################################
# Java
#############################################################################

# If you don't need lcio.jar you can set this to OFF
#SET( INSTALL_JAR OFF CACHE BOOL "Set to OFF to skip build/install lcio.jar" FORCE )


#############################################################################
# Project options
#############################################################################

#SET( INSTALL_DOC OFF CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )
#SET( BUILD_LCIO_TESTJOBS ON CACHE BOOL "Set to ON to build LCIO testjobs" FORCE )
#SET( BUILD_F77_TESTJOBS ON CACHE BOOL "Set to ON to build LCIO F77 testjobs" FORCE )

# set cmake build type, default value is: RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
#SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE )

#############################################################################
# Advanced options
#############################################################################

#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

# installation path for LCIO
#SET( CMAKE_INSTALL_PREFIX "/foo/bar" CACHE STRING "Where to install LCIO" FORCE )

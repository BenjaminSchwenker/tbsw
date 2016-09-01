#############################################################################
# cmake build setup for LCCD
#
# For building lccd with cmake type:
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


# Path to LCIO
SET( LCIO_HOME "$ENV{TBSW_HOME}/lcio"
    CACHE PATH "Path to LCIO" FORCE )

# CMake Modules Path
SET( CMAKE_MODULE_PATH "$ENV{TBSW_HOME}/CMakeModules"
    CACHE PATH "Path to CMake Modules" FORCE )





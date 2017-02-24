#############################################################################
#
# CMAKE build setup for EudaqInput
#
# For building Marlin with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
#
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

# CMake Modules Path
SET( CMAKE_MODULE_PATH "$ENV{TBSW_HOME}/CMakeModules"
    CACHE PATH "Path to CMake Modules" FORCE )

# Path to LCIO
SET( LCIO_HOME "$ENV{TBSW_HOME}/lcio"
    CACHE PATH "Path to LCIO" FORCE )

# Path to Marlin
SET( Marlin_HOME "$ENV{TBSW_HOME}/Marlin"
     CACHE PATH "Path to Marlin" FORCE )


#############################################################################
# Project Options
#############################################################################

# Project depends on ...
SET( PROJECT_DEPENDS "LCIO Marlin" 
     CACHE STRING "EudaqInput dependence" FORCE )
    
# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel), default: RelWithDebInfo
SET( CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build" FORCE )

# Set basic options
SET( INSTALL_DOC OFF CACHE BOOL "Set to ON not to skip build/install Documentation" FORCE )




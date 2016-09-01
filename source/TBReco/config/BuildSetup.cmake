#############################################################################
#
# CMAKE build setup for TBReco
#
# For building TBReco with cmake type:
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

# Path to GEAR
SET( GEAR_HOME "$ENV{TBSW_HOME}/gear"
    CACHE PATH "Path to GEAR" FORCE )

# Path to Marlin
SET( Marlin_HOME "$ENV{TBSW_HOME}/Marlin"
     CACHE PATH "Path to Marlin" FORCE )

# Path to TBTOOLS
SET( TBTools_HOME "$ENV{TBSW_HOME}/TBTools"
     CACHE PATH "Path to TBTools" FORCE )

# Path to CLHEP, respectively to 'clhep-config'
SET( CLHEP_HOME "$ENV{CLHEP_HOME}"  
     CACHE PATH "Path to CLHEP" FORCE )

# Path to ROOT
SET( ROOT_HOME "$ENV{ROOTSYS}"
    CACHE PATH "Path to ROOT " FORCE )

#############################################################################
# Project Options
#############################################################################

# Project depends on ...
SET( PROJECT_DEPENDS "CLHEP GEAR LCIO Marlin ROOT TBTools" 
     CACHE STRING "TBReco dependence" FORCE )
    
# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel), default: RelWithDebInfo
SET( CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build" FORCE )

# Set basic options
SET( INSTALL_DOC OFF CACHE BOOL "Set to ON not to skip build/install Documentation" FORCE )




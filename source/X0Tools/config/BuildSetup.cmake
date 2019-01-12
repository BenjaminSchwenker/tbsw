#############################################################################
#
# CMAKE build setup for X0Tools
#
# For building X0Tools with cmake type:
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

# Path to ROOT
SET( ROOT_HOME "$ENV{ROOTSYS}"
    CACHE PATH "Path to ROOT " FORCE )

#############################################################################
# Project Options
#############################################################################

# Project depends on ...
SET( PROJECT_DEPENDS "ROOT" 
     CACHE STRING "X0Tools dependence" FORCE )
    
# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel), default: RelWithDebInfo
SET( CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build" FORCE )

# Set basic options
SET( INSTALL_DOC OFF CACHE BOOL "Set to ON not to skip build/install Documentation" FORCE )




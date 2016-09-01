#############################################################################
#
# CMAKE build setup for GEAR
#
# For building Gear with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
#
#############################################################################


#############################################################################
# Setup project options
#############################################################################

SET( BUILD_32BIT_COMPATIBLE OFF 
    CACHE BOOL "Set to Off for 64bit build" FORCE   )


# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel), default: RelWithDebInfo
SET( CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build" FORCE )

# Set basic options
SET( BUILD_GEAR_TESTS ON CACHE BOOL  "Set to OFF not to build GEAR tests" FORCE )
SET( INSTALL_DOC ON CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )

# Set advanced options
#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )


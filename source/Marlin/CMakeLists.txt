########################################################
# CMake compatibility issues: don't modify this, please!
#CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( Marlin )

# project version
SET( Marlin_MAJOR_VERSION 0 )
SET( Marlin_MINOR_VERSION 11 )
SET( Marlin_PATCH_LEVEL 0 )

# definitions
IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
    CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
    ADD_DEFINITIONS( "-DMARLINDEBUG" )
ENDIF()
                                                                                                                                                           
IF( MARLIN_NO_DLL )
    ADD_DEFINITIONS( "-DMARLIN_NO_DLL" )
ENDIF()

##########################################################################################

# library *nix style versioning
SET( Marlin_SOVERSION "${Marlin_MAJOR_VERSION}.${Marlin_MINOR_VERSION}" )
SET( Marlin_VERSION "${Marlin_SOVERSION}.${Marlin_PATCH_LEVEL}" )



# add install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# add install path to the rpath list (apple)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )



# input directories
ADD_SUBDIRECTORY( streamlog )
ADD_SUBDIRECTORY( src )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------- ${PROJECT_NAME} --------------------------------------" )
MESSAGE( STATUS "MARLIN_NO_DLL = ${MARLIN_NO_DLL}" )
MESSAGE( STATUS "Marlin_VERSION = ${Marlin_VERSION}" )


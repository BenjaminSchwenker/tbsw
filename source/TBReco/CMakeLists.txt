#############################################################################
#
# CMAKE file for building TBReco
# 
# @author Benjamin Schwenker
#
#############################################################################

# Project name
PROJECT( TBReco )

# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


### SETTINGS ################################################################

# Project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )



### ADVANCED SETTINGS #######################################################

# Library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )


### INCLUDE, SOURCE, LIB ####################################################

# Include directories
INCLUDE_DIRECTORIES( include ${Eigen_INCLUDE_PATH} ${ROOT_INCLUDE_DIRS})

# Set list of sources
AUX_SOURCE_DIRECTORY( src lib_srcs)
file(GLOB_RECURSE lib_hdrs RELATIVE "${PROJECT_SOURCE_DIR}" *.h *.hh)

# Create library with the list of library sources
ADD_LIBRARY( ${PROJECT_NAME} ${lib_srcs} ${lib_hdrs} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} lib_LCIO lib_Marlin TBTools)

# Install libraries
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

### FINAL STATUS ############################################################

# Display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------- ${PROJECT_NAME} --------------------------------------" )
MESSAGE( STATUS "${PROJECT_NAME}_VERSION = ${${PROJECT_NAME}_VERSION}" )





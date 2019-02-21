#############################################################################
#
# CMAKE file for building EudaqInput
# 
# @author Benjamin Schwenker
#
#############################################################################


# CMake compatibility issues: don't modify this, please!
#CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

### PROJECT #################################################################

# Project name
PROJECT( EudaqInput )

# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


### SETTINGS ################################################################

# Project options
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )



### ADVANCED SETTINGS #######################################################

# Library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )


### INCLUDE, SOURCE, LIB ####################################################
# Include directories
INCLUDE_DIRECTORIES( include )
file(GLOB_RECURSE lib_hdrs RELATIVE "${PROJECT_SOURCE_DIR}" *.h *.hh)

# Set list of sources
AUX_SOURCE_DIRECTORY( src lib_srcs )

# Create library with the list of library sources
ADD_LIBRARY( ${PROJECT_NAME} ${lib_srcs} ${lib_hdrs})
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} lib_LCIO lib_Marlin)

# Install libraries
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )
# display status message for important variables
if (VERBOSE)
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------- ${PROJECT_NAME} --------------------------------------" )
MESSAGE( STATUS "${PROJECT_NAME}_VERSION = ${${PROJECT_NAME}_VERSION}" )
ENDIF()

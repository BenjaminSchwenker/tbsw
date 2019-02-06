#################################
# cmake file for building GEAR
# @author Jan Engels, DESY
# @author Benjamin Schwenker, Uni Göttingen
#################################

########################################################
# CMake compatibility issues: don't modify this, please!
#CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( GEAR )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 11 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 1 )

# library *nix style versioning
SET( GEAR_SOVERSION "${GEAR_MAJOR_VERSION}.${GEAR_MINOR_VERSION}" )
SET( GEAR_VERSION "${GEAR_SOVERSION}.${GEAR_PATCH_LEVEL}" )

# create symbolic link for include files
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/src/include/gearxml" )
    MESSAGE( STATUS "Creating symlink: ${PROJECT_SOURCE_DIR}/src/include/gearxml..." )
    EXEC_PROGRAM( "${CMAKE_COMMAND}"
        "${PROJECT_SOURCE_DIR}/src/include"
        ARGS -E create_symlink "${PROJECT_SOURCE_DIR}/src/gearxml/include" "gearxml"
        OUTPUT_VARIABLE out_tmp )
ENDIF()
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/src/include/gearimpl" )
    MESSAGE( STATUS "Creating symlink: ${PROJECT_SOURCE_DIR}/src/include/gearimpl..." )
    EXEC_PROGRAM( "${CMAKE_COMMAND}"
        "${PROJECT_SOURCE_DIR}/src/include"
        ARGS -E create_symlink "${PROJECT_SOURCE_DIR}/src/gearimpl/include" "gearimpl"
        OUTPUT_VARIABLE out_tmp )
ENDIF()

##########################################################################################


# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

SET(gear_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/include")

INCLUDE_DIRECTORIES( "src/include" )

# input directories
ADD_SUBDIRECTORY( src/gearxml )
ADD_SUBDIRECTORY( src/gearimpl )






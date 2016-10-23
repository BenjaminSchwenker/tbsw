#############################################################################
#
# CMAKE file for building TBTools
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
PROJECT( TBTools )

# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


### SETTINGS ################################################################

# Project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

# Project dependencies
IF( PROJECT_DEPENDS )
    SET( ${PROJECT_NAME}_DEPENDS "${PROJECT_DEPENDS}" )
ELSE()
    SET( ${PROJECT_NAME}_DEPENDS "CLHEP GEAR LCIO Marlin ROOT")
ENDIF()    

# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel)
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# Set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" )
ENDIF()

#---------------- 32/64 bit issues ---------------------------------------
IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "64 bit architecture detected" )

    IF( BUILD_32BIT_COMPATIBLE )
        IF( COMMAND SET_PROPERTY )
            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
            SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
        ELSE()
            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
        ENDIF()
        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
    ELSE()
        MESSAGE( STATUS "Build in native 64 bit mode" )
    ENDIF()
ENDIF()
#-------------------------------------------------------------------------


### ADVANCED SETTINGS #######################################################

# Library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

# Output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# Add install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# Add install path to the rpath list (APPLE)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# Append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )


### DEPENDENCIES #######################################################

# Use Module Macro to find out if all dependencies are Ok
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    # load macro
    IF( NOT EXISTS "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
        MESSAGE( FATAL_ERROR
	    "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake " 
            "\nSorry, could not find MacroCheckDeps.cmake...\n"
            "Please set CMAKE_MODULE_PATH correctly with: "
            "cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
    ENDIF()
    INCLUDE( "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
    CHECK_DEPS()
ENDIF()


### create root dictionary #################################################

# create a new dictionary for class AlignEvent
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/source/include/AlignEventDict.h" )
    MESSAGE( STATUS "Creating root dictionary for class AlignEvent" )
    
    EXEC_PROGRAM( rootcint "${PROJECT_SOURCE_DIR}/source/include"
            ARGS -f AlignEventDict.C 
            ARGS -c AlignEvent.h LinkDef.h
            OUTPUT_VARIABLE out_tmp
            RETURN_VALUE out_ret )
    IF( out_ret )
        MESSAGE( STATUS "Error executing rootcint!!" )
    ENDIF()
    
ENDIF()

# move dictionary file to other sources 
IF( EXISTS "${PROJECT_SOURCE_DIR}/source/include/AlignEventDict.C" )
    MESSAGE( STATUS "Created dictionary for class AlignEvent" )
    
    EXEC_PROGRAM( mv "${PROJECT_SOURCE_DIR}/source/include"
            ARGS AlignEventDict.C 
            ARGS ${PROJECT_SOURCE_DIR}/source/src
            OUTPUT_VARIABLE out_tmp
            RETURN_VALUE out_ret )
    IF( out_ret )
        MESSAGE( STATUS "Error moving dictionary!!" )
    ENDIF()
    
ENDIF()

### INCLUDE, SOURCE, LIB ####################################################

# Compiler flags
#ADD_DEFINITIONS( "-Wall -Wno-long-long -ansi -pedantic -std=c++11" )
ADD_DEFINITIONS( "-Wall -Wno-long-long -ansi -pedantic -std=c++1y" )


# Include directories
INCLUDE_DIRECTORIES( source/include )

# Set list of sources
AUX_SOURCE_DIRECTORY( source/src lib_srcs)

IF( COMMAND cmake_policy )
  cmake_policy(SET CMP0003 NEW)
ENDIF( COMMAND cmake_policy )


# Create library with the list of library sources
ADD_LIBRARY( ${PROJECT_NAME} ${lib_srcs} )
#ADD_LIBRARY( ${PROJECT_NAME} source/src/RombIntSolver.cc )

# Create target for libraries
ADD_CUSTOM_TARGET( lib DEPENDS ${PROJECT_NAME} )

# Install libraries
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

# Create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/config/${PROJECT_NAME}.cmake.in"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY )

# Install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION . )


# Create uninstall configuration file 
CONFIGURE_FILE( "${CMAKE_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
                "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# Create uninstall target
ADD_CUSTOM_TARGET( uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake" )

### FINAL STATUS ############################################################

# Display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
MESSAGE( STATUS "${PROJECT_NAME}_DEPENDS = \"${${PROJECT_NAME}_DEPENDS}\"" )
MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )


# Force some variables to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install ${PROJECT_NAME}" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )


# Export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "${PROJECT_NAME}BuildSettings.cmake" )

# Export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )


#################################
# cmake file for building LCCD
# @author Jan Engels, DESY
#################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( LCCD )

# project version
SET( LCCD_MAJOR_VERSION 0 )
SET( LCCD_MINOR_VERSION 6 )
SET( LCCD_PATCH_LEVEL 1 )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( BUILD_LCCD_TESTS "Set to ON to build LCCD tests" OFF )
OPTION( BUILD_CONDDB_TESTS "Set to ON to build CondDB tests" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

# project dependencies
SET( LCCD_DEPENDS "LCIO" )

# set default cmake build type to RelWithDebInfo
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${LCCD_SOURCE_DIR}" )
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

# specific LCCD variable
IF( NOT LCCD_DB_INIT )
    SET( LCCD_DB_INIT "localhost:lccd_test:calvin:hobbes" )
ENDIF()

# omit warning about preprocessor definitions being escaped
IF(COMMAND CMAKE_POLICY)
    CMAKE_POLICY(SET CMP0005 OLD)
ENDIF()

#----- need long long for int64 for now ------
#FIXME: should establish wether we are on a 32bit or 64 bit machine ....
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-DLCCD_DB_INIT_DEFAULT=\\\"${LCCD_DB_INIT}\\\"" )

# library *nix style versioning
SET( LCCD_SOVERSION "${LCCD_MAJOR_VERSION}.${LCCD_MINOR_VERSION}" )
SET( LCCD_VERSION "${LCCD_SOVERSION}.${LCCD_PATCH_LEVEL}" )

##########################################################################################

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

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${LCCD_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${LCCD_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# DEPENDENCIES: this code has to be placed before adding any library or
# executable so that these are linked properly against the dependencies
IF( DEFINED LCCD_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    IF( NOT DEFINED CMAKE_MODULE_PATH )
        MESSAGE( FATAL_ERROR "CMAKE_MODULE_PATH not set! Set it with: "
            "-DCMAKE_MODULE_PATH=\"/path_to_ilcsoft/CMakeModules\"" )
    ENDIF()
    # load macro
    INCLUDE( "MacroCheckDeps" )
    CHECK_DEPS()
ENDIF()

# input directories
ADD_SUBDIRECTORY( source )

# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${LCCD_SOURCE_DIR}/LCCDConfig.cmake.in"
                "${LCCD_BINARY_DIR}/LCCDConfig.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${LCCD_BINARY_DIR}/LCCDConfig.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
ENDIF()
MESSAGE( STATUS "BUILD_LCCD_TESTS = ${BUILD_LCCD_TESTS}" )
MESSAGE( STATUS "BUILD_CONDDB_TESTS = ${BUILD_CONDDB_TESTS}" )
MESSAGE( STATUS "LCCD_DB_INIT = ${LCCD_DB_INIT}" )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install LCCD" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )
SET( BUILD_LCCD_TESTS "${BUILD_LCCD_TESTS}" CACHE BOOL
    "Set to ON to build LCCD tests" FORCE )
SET( BUILD_CONDDB_TESTS "${BUILD_CONDDB_TESTS}" CACHE BOOL
    "Set to ON to build CONDDB tests" FORCE )
SET( LCCD_DB_INIT "${LCCD_DB_INIT}" CACHE STRING
    "DB initialization" FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "LCCDBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/LCCDBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "LCCDLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/LCCDLibDeps.cmake" DESTINATION lib/cmake )

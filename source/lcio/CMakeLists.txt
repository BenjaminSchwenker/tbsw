#################################
# cmake file for building LCIO
# @author Frank Gaede, DESY
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
# enable testing
ENABLE_TESTING()
########################################################

# project name
PROJECT( LCIO )

# project version
SET( LCIO_MAJOR_VERSION 1 )
SET( LCIO_MINOR_VERSION 12 )
SET( LCIO_PATCH_LEVEL 2 )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( BUILD_LCIO_TESTS "Set to ON to build LCIO tests under src/cpp/src/TESTS" OFF )
OPTION( BUILD_LCIO_EXAMPLES "Set to ON to build LCIO examples under src/cpp/src/EXAMPLES" OFF )
OPTION( BUILD_WITH_DCAP "Set to ON to build LCIO with dcap (dCache) support" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )
OPTION( INSTALL_JAR "Set to OFF to skip build/install lcio.jar" OFF )
OPTION( BUILD_ROOTDICT "Set to ON to build ROOT dictionary" OFF )


# set default cmake build type to RelWithDebInfo
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

#---------------- 32/64 bit issues ---------------------------------------
#INCLUDE( CheckTypeSize )
#CHECK_TYPE_SIZE( "void*" SIZE_OF_VOID_P )

IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "64 bit architecture detected" )

    IF( BUILD_32BIT_COMPATIBLE )
        IF( COMMAND SET_PROPERTY )
            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
            SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
            SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32" )
            SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32" )
            SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -m32" )
            SET( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -m32" )
        ELSE()
            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
        ENDIF()
        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
    ELSE()
        MESSAGE( STATUS "Build in native 64 bit mode" )
    ENDIF()
ENDIF()



##########################################################################################

# library *nix style versioning
SET( LCIO_SOVERSION "${LCIO_MAJOR_VERSION}.${LCIO_MINOR_VERSION}" )
SET( LCIO_VERSION "${LCIO_SOVERSION}.${LCIO_PATCH_LEVEL}" )

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

# build output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# input directories
IF( BUILD_WITH_DCAP )
  IF( NOT DEFINED DCAP_HOME )
      MESSAGE( FATAL_ERROR "Check for DCAP: DCAP_HOME not set!! "
              "Please call cmake with -DDCAP_HOME=<path_to_DCAP>" )
  ENDIF()
  INCLUDE_DIRECTORIES( ${DCAP_HOME}/include )
  ADD_DEFINITIONS( "-DSIO_USE_DCAP" )
  FIND_LIBRARY( DCAP_LIB NAMES dcap PATHS ${DCAP_HOME}/lib NO_DEFAULT_PATH )
  LINK_LIBRARIES( ${DCAP_LIB} ) 
ENDIF()

ADD_SUBDIRECTORY( sio sio )
ADD_SUBDIRECTORY( src/cpp lcio )


# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/LCIOConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/LCIOConfig.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/LCIOConfig.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
ENDIF()
MESSAGE( STATUS "BUILD_LCIO_TESTS = ${BUILD_LCIO_TESTS}" )
MESSAGE( STATUS "BUILD_LCIO_EXAMPLES = ${BUILD_LCIO_EXAMPLES}" )
MESSAGE( STATUS "BUILD_WITH_DCAP = ${BUILD_WITH_DCAP}" )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "INSTALL_JAR = ${INSTALL_JAR}" )
MESSAGE( STATUS "BUILD_ROOTDICT = ${BUILD_ROOTDICT}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install LCIO" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( BUILD_LCIO_TESTS "${BUILD_LCIO_TESTS}" CACHE BOOL
    "Set to ON to build LCIO tests under src/cpp/src/TESTS" FORCE )
SET( BUILD_LCIO_EXAMPLES "${BUILD_LCIO_EXAMPLES}" CACHE BOOL
    "Set to ON to build LCIO examples under src/cpp/src/EXAMPLES" FORCE )
SET( BUILD_WITH_DCAP "${BUILD_WITH_DCAP}" CACHE BOOL
    "Set to ON to build LCIO with dcap support (dCache)" FORCE )
SET( DCAP_HOME "${DCAP_HOME}" CACHE PATH
    "Set to dcap installation directory - only if BUILD_WITH_DCAP=ON" FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )
SET( INSTALL_JAR "${INSTALL_JAR}" CACHE BOOL
    "Set to OFF to skip build/install lcio.jar" FORCE )
SET( BUILD_ROOTDICT "${BUILD_ROOTDICT}" CACHE BOOL
     "Set to ON to build LCIO ROOT dictionary" FORCE )


# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "LCIOBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/LCIOBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "LCIOLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/LCIOLibDeps.cmake" DESTINATION lib/cmake )

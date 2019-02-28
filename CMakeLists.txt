cmake_minimum_required(VERSION 2.4.6)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( TBSW )

# project version
SET( tbsw_MAJOR_VERSION 0 )
SET( tbsw_MINOR_VERSION 11 )
SET( tbsw_PATCH_LEVEL 0 )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )
OPTION( BUILD_EudaqInput "Set to OFF to not build EudaqInput" ON )
OPTION( BUILD_TBTools "Set to OFF to not build TBTools" ON )
OPTION( BUILD_TBReco "Set to OFF to not build TBReco" ON )
OPTION( BUILD_X0Tools "Set to OFF to build X0Tools" ON )
OPTION( VERBOSE "Set to ON to make Verbose makefile" OFF )
# LCIO Options
OPTION( BUILD_LCIO_TESTS "Set to ON to build LCIO tests under src/cpp/src/TESTS" OFF )
OPTION( BUILD_LCIO_EXAMPLES "Set to ON to build LCIO examples under src/cpp/src/EXAMPLES" OFF )
OPTION( BUILD_WITH_DCAP "Set to ON to build LCIO with dcap (dCache) support" OFF )
OPTION( INSTALL_JAR "Set to OFF to skip build/install lcio.jar" OFF )
OPTION( BUILD_ROOTDICT "Set to ON to build ROOT dictionary" OFF )
# MARLIN Options
OPTION( MARLIN_NO_DLL "Set to ON to build Marlin without DLL support" OFF )


# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel)
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

IF(ROOT_HOME)
    SET( ROOT_BASH_PATH ${ROOT_HOME})
    SET( CMAKE_PREFIX_PATH ${ROOT_HOME})
ELSE()
    SET( ROOT_BASH_PATH $ENV{ROOTSYS})
    SET( CMAKE_PREFIX_PATH $ENV{ROOTSYS})
ENDIF()
SET( CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE PATH
    "CMAKE_PREFIX_PATH" FORCE )
find_package(ROOT REQUIRED)
# display status message for important variables
IF( VERBOSE )
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------------- ROOT -------------------------------------" )
MESSAGE( STATUS "ROOT_FOUND        = ${ROOT_FOUND}" )
MESSAGE( STATUS "ROOT_INCLUDE_DIRS  = ${ROOT_INCLUDE_DIRS}" )
MESSAGE( STATUS "ROOT_rootcint_CMD  = ${ROOT_rootcint_CMD}" )
MESSAGE( STATUS "ROOT_LIBRARIES  = ${ROOT_LIBRARIES}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH  = ${CMAKE_PREFIX_PATH}" )
ENDIF()




#-------------Compiler Flags---------------------------------
#------ all compiler flags should go here ----
ADD_DEFINITIONS( "-std=c++14" )
#------ correct C++
ADD_DEFINITIONS( "-Wall -pedantic" )
#------ performance
ADD_DEFINITIONS( "-march=native -mtune=native" )
#ADD_DEFINITIONS( "-mtune=native -mno-avx -mno-avx2" )
#------ lets get real pedantic as a test...
ADD_DEFINITIONS( "-Wextra -Wno-unused-macros -Wno-shadow ")# -Wdouble-promotion ")
#------ but not too much...
ADD_DEFINITIONS( "-Wno-unused-function ")
# definitions
IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
    CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
    ADD_DEFINITIONS( "-g" )
ENDIF()
IF( VERBOSE )
    set(CMAKE_VERBOSE_MAKEFILE ON)
ENDIF()

if(POLICY CMP0002)
  cmake_policy(SET CMP0002 NEW)
endif()
IF(POLICY CMP0003)
  cmake_policy(SET CMP0003 NEW)
ENDIF()
IF(POLICY CMP0009)
  cmake_policy(SET CMP0009 NEW)
ENDIF()

# set default install prefix to project root directory

IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}" )
ENDIF()

#---------------- 32/64 bit issues ---------------------------------------
IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    IF( BUILD_32BIT_COMPATIBLE )
        IF( COMMAND SET_PROPERTY )
            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
            SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
        ELSE()
            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
        ENDIF()
        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
    ELSE()
        SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64" )
        MESSAGE( STATUS "Build in native 64 bit mode" )
    ENDIF()
ENDIF()

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

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



IF(NOT Eigen_INCLUDE_PATH )
    SET(  Eigen_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/source/Eigen3Include")
ENDIF()

# input directories
ADD_SUBDIRECTORY( source )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "---------------------------- TBSW Summary-----------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS      = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE       = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "ROOT_INCLUDE_DIRS      = ${ROOT_INCLUDE_DIRS}" )
MESSAGE( STATUS "VERBOSE                = ${VERBOSE}" )
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
ENDIF()
MESSAGE( STATUS "BUILD_EudaqInput       = ${BUILD_EudaqInput}" )
MESSAGE( STATUS "BUILD_TBTools          = ${BUILD_TBTools}" )
MESSAGE( STATUS "BUILD_TBReco           = ${BUILD_TBReco}" )
MESSAGE( STATUS "BUILD_X0Tools          = ${BUILD_X0Tools}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install Marlin" FORCE )
SET( Eigen_INCLUDE_PATH "${Eigen_INCLUDE_PATH}" CACHE PATH
    "Eigen Library header files" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( MARLIN_NO_DLL "${MARLIN_NO_DLL}" CACHE BOOL
    "Set to ON to build Marlin without DLL support" FORCE )

SET( BUILD_EudaqInput "${BUILD_EudaqInput}" CACHE BOOL
     "Set to OFF to not build EudaqInput"  FORCE )
SET( BUILD_TBTools "${BUILD_TBTools}" CACHE BOOL
     "Set to OFF to not build TBTools"  FORCE )
SET( BUILD_TBReco "${BUILD_TBReco}" CACHE BOOL
     "Set to OFF to not build TBReco"  FORCE )
SET( BUILD_X0Tools "${BUILD_X0Tools}" CACHE BOOL
      "Set to ON to build TBReco"  FORCE )

# Create the init_tbsw.sh file

file(WRITE  "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#!/bin/bash\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#--------------------------------------------------------------------------------\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#    ROOT                                                                        \n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#--------------------------------------------------------------------------------\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export ROOTSYS=${ROOT_BASH_PATH}\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export PATH=${ROOT_BASH_PATH}/bin:\$PATH\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export LD_LIBRARY_PATH=${ROOT_BASH_PATH}/lib:\$LD_LIBRARY_PATH\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export PYTHONPATH=${ROOT_BASH_PATH}/lib:${PROJECT_SOURCE_DIR}/source:\$PYTHONPATH\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export ROOT_INCLUDE_PATH=${PROJECT_BINARY_DIR}\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "\n")

file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#--------------------------------------------------------------------------------\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#    TBSW                                                                        \n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "#--------------------------------------------------------------------------------\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export PATH=${PROJECT_BINARY_DIR}/bin:\$PATH\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/lib:\$LD_LIBRARY_PATH\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export MARLIN_DLL=${PROJECT_BINARY_DIR}/lib/libTBReco.so:${PROJECT_BINARY_DIR}/lib/libEudaqInput.so:\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "export MARLIN=${PROJECT_BINARY_DIR}\n")
file(APPEND "${PROJECT_BINARY_DIR}/init_tbsw.sh" "\n")

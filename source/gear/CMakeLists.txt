#################################
# cmake file for building GEAR
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
PROJECT( GEAR )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 11 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 1 )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

# set default cmake build type to RelWithDebInfo
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
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



# library *nix style versioning
SET( GEAR_SOVERSION "${GEAR_MAJOR_VERSION}.${GEAR_MINOR_VERSION}" )
SET( GEAR_VERSION "${GEAR_SOVERSION}.${GEAR_PATCH_LEVEL}" )

# create symbolic link for include files
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/src/cpp/include/gearxml" )
    MESSAGE( STATUS "Creating symlink: ${PROJECT_SOURCE_DIR}/src/cpp/include/gearxml..." )
    EXEC_PROGRAM( "${CMAKE_COMMAND}"
        "${PROJECT_SOURCE_DIR}/src/cpp/include"
        ARGS -E create_symlink "${PROJECT_SOURCE_DIR}/src/cpp/gearxml/include" "gearxml"
        OUTPUT_VARIABLE out_tmp )
ENDIF()
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/src/cpp/include/gearimpl" )
    MESSAGE( STATUS "Creating symlink: ${PROJECT_SOURCE_DIR}/src/cpp/include/gearimpl..." )
    EXEC_PROGRAM( "${CMAKE_COMMAND}"
        "${PROJECT_SOURCE_DIR}/src/cpp/include"
        ARGS -E create_symlink "${PROJECT_SOURCE_DIR}/src/cpp/gearimpl/include" "gearimpl"
        OUTPUT_VARIABLE out_tmp )
ENDIF()

# generate header files with ant aid.generate
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/src/cpp/include/gear/GEAR.h" )
    MESSAGE( STATUS "Generating GEAR header files with ant aid.generate..." )
    MESSAGE( STATUS "Autodetecting java..." )
    FIND_PACKAGE( Java )
    IF( NOT JAVA_RUNTIME OR NOT JAVA_COMPILE )
        MESSAGE( FATAL_ERROR "Java not found!! Please make sure that java is in your PATH." )
    ENDIF()
    # extract java bin path out of full path to java runtime
    STRING( REGEX REPLACE "(.*)\\/java$" "\\1" JAVA_RUNTIME_PATH ${JAVA_RUNTIME} )
    # parse the output of java -version
    EXEC_PROGRAM( "${JAVA_RUNTIME}" ARGS "-version"
            OUTPUT_VARIABLE out_tmp
            RETURN_VALUE out_ret )
    IF( out_ret )
        MESSAGE( STATUS "Error executing java -version!!" )
    ENDIF()
    STRING( REGEX REPLACE ".* version \"([0-9]+\\.[0-9]+\\.[0-9]+).*"
            "\\1" JAVA_VERSION "${out_tmp}" )

    # display info
    MESSAGE( STATUS "Detected Java version ${JAVA_VERSION}..." )

    # set environment variables for
    # generating header files with ant aid.generate
    SET( ENV{GEAR} "${PROJECT_SOURCE_DIR}" )
    SET( ENV{PATH} "${PROJECT_SOURCE_DIR}/tools:${JAVA_RUNTIME_PATH}:$ENV{PATH}" )

    # generate header files with ant aid.generate
    EXEC_PROGRAM( ant "${PROJECT_SOURCE_DIR}"
        ARGS aid.generate
        OUTPUT_VARIABLE out_tmp
        RETURN_VALUE out_ret )
    # check if everything went ok
    IF( NOT out_ret )
        MESSAGE( STATUS "Header files generated in ${PROJECT_SOURCE_DIR}/src/cpp/include/gear..." )
    ELSE()
        MESSAGE( STATUS "Error generating header files!!" )
        MESSAGE( FATAL_ERROR "${out_tmp}" )
    ENDIF()
ENDIF()


# code for *nix only!
IF( UNIX )
    # find shell
    FIND_PROGRAM( SH
      sh
      ${CYGWIN_INSTALL_PATH}/bin
      /bin
      /usr/bin
      /usr/local/bin
      /sbin
      # FIXME add path to MacOS here
    )
    MARK_AS_ADVANCED( SH )

    IF( SH )
        # generate shell script for setting environment before running a command
        MESSAGE( STATUS "Generating shell environment script..." )
        FILE( WRITE "${PROJECT_BINARY_DIR}/script.sh"
                    "export GEAR=${PROJECT_SOURCE_DIR}\n"
                    "export GEARVERSION=${GEAR_VERSION}\n"
                    "export PATH=${JAVA_RUNTIME_PATH}/bin:${PROJECT_SOURCE_DIR}/tools:$PATH\n"
                    "$1 $2" )

        # create "always out-of-date" headers target
        ADD_CUSTOM_TARGET( headers
            COMMAND ${SH} "${PROJECT_BINARY_DIR}/script.sh" ant aid.generate
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            COMMENT "Generating header files with ant aid.generate..."
            VERBATIM )
        
    ELSE()
        MESSAGE( STATUS "Unix Shell not found - cannot enable headers target!!" )
        IF( INSTALL_DOC )
            MESSAGE( STATUS "Unix Shell not found - INSTALL_DOC forced to OFF" )
            SET( INSTALL_DOC OFF )
        ENDIF()
    ENDIF()
ENDIF()



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
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib )

# FIXME: check for 32/64 bit
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

INCLUDE_DIRECTORIES( "src/cpp/include" )

# input directories
ADD_SUBDIRECTORY( src/cpp/gearxml )
ADD_SUBDIRECTORY( src/cpp/gearimpl )
ADD_SUBDIRECTORY( src/cpp/test )

# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/GEARConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/GEARConfig.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/GEARConfig.cmake" DESTINATION . )

# install include files
INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/src/cpp/include"
        DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
ENDIF()
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "To change a value type: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install GEAR" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "GEARBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/GEARBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "GEARLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/GEARLibDeps.cmake" DESTINATION lib/cmake )


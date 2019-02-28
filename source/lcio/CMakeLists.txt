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


# set default cmake build type to RelWithDebInfo
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()


##########################################################################################

# library *nix style versioning
SET( LCIO_SOVERSION "${LCIO_MAJOR_VERSION}.${LCIO_MINOR_VERSION}" )
SET( LCIO_VERSION "${LCIO_SOVERSION}.${LCIO_PATCH_LEVEL}" )

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
ADD_SUBDIRECTORY( lcio lcio )


# display status message for important variables
if (VERBOSE)
MESSAGE( STATUS )
MESSAGE( STATUS "----------------------------- LCIO --------------------------------------------" )
MESSAGE( STATUS "BUILD_LCIO_TESTS = ${BUILD_LCIO_TESTS}" )
MESSAGE( STATUS "BUILD_LCIO_EXAMPLES = ${BUILD_LCIO_EXAMPLES}" )
MESSAGE( STATUS "BUILD_WITH_DCAP = ${BUILD_WITH_DCAP}" )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "INSTALL_JAR = ${INSTALL_JAR}" )
MESSAGE( STATUS "BUILD_ROOTDICT = ${BUILD_ROOTDICT}" )
ENDIF()

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_LCIO_TESTS "${BUILD_LCIO_TESTS}" CACHE BOOL
    "Set to ON to build LCIO tests under src/cpp/src/TESTS" FORCE )
SET( BUILD_LCIO_EXAMPLES "${BUILD_LCIO_EXAMPLES}" CACHE BOOL
    "Set to ON to build LCIO examples under src/cpp/src/EXAMPLES" FORCE )
SET( BUILD_WITH_DCAP "${BUILD_WITH_DCAP}" CACHE BOOL
    "Set to ON to build LCIO with dcap support (dCache)" FORCE )
SET( DCAP_HOME "${DCAP_HOME}" CACHE PATH
    "Set to dcap installation directory - only if BUILD_WITH_DCAP=ON" FORCE )
SET( INSTALL_JAR "${INSTALL_JAR}" CACHE BOOL
    "Set to OFF to skip build/install lcio.jar" FORCE )
SET( BUILD_ROOTDICT "${BUILD_ROOTDICT}" CACHE BOOL
     "Set to ON to build LCIO ROOT dictionary" FORCE )



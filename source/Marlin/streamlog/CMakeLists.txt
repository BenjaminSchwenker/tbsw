#################################
# cmake file for building streamlog
# @author Frank Gaede, DESY
#################################

# cmake minimum required version
#CMAKE_MINIMUM_REQUIRED( VERSION 2.4 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.2 )
#MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# project name
PROJECT( streamlog )

# project options
OPTION( BUILD_SHARED_LIBS "Set this to OFF to build static libraries " ON )
#OPTION( BUILD_streamlog_TESTS "Set to ON to build streamlog tests" OFF )

# project version
SET( streamlog_MAJOR_VERSION 0 )
SET( streamlog_MINOR_VERSION 1 )
SET( streamlog_PATCH_LEVEL 0 )

SET( streamlog_SOVERSION "${streamlog_MAJOR_VERSION}.${streamlog_MINOR_VERSION}" )
SET( streamlog_VERSION "${streamlog_SOVERSION}.${streamlog_PATCH_LEVEL}" )

# include directories
INCLUDE_DIRECTORIES( BEFORE ./source/include )

# get list of all source files
AUX_SOURCE_DIRECTORY( ./source/src lib_srcs )

# create library with the list of library sources
ADD_LIBRARY( lib_streamlog ${lib_srcs} ${streamlog_hdrs} )

# create symbolic lib target for calling target lib_XXX
#ADD_CUSTOM_TARGET( lib DEPENDS lib_streamlog )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_streamlog PROPERTIES
	# create *nix style library versions + symbolic links
	VERSION ${streamlog_VERSION}
	SOVERSION ${streamlog_SOVERSION}
	# allow building static and shared libraries on the same build tree
	CLEAN_DIRECT_OUTPUT 1
	# avoid conflicts between library and binary target names
	OUTPUT_NAME streamlog)

# install library
INSTALL( TARGETS lib_streamlog DESTINATION lib PERMISSIONS
		OWNER_READ OWNER_WRITE OWNER_EXECUTE
		GROUP_READ GROUP_EXECUTE
		WORLD_READ WORLD_EXECUTE )

#-----------------------------------------------------------
ADD_SUBDIRECTORY( examples test EXCLUDE_FROM_ALL )


# install include files
INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/source/include/streamlog"
    DESTINATION include
    PATTERN "*~" EXCLUDE
    PATTERN "*CVS*" EXCLUDE )



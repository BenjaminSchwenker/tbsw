#################################
# cmake file for building Marlin
# @author Jan Engels, DESY
# @author Benjamin Schwenker, Uni Goettingen, removed MarlinGUI because not needed in tbsw
#################################

####################################################################################################
### INCLUDES
####################################################################################################

SET( Marlin_Version_Header "MarlinVersion.h" )

FILE( WRITE "${PROJECT_BINARY_DIR}/marlin/${Marlin_Version_Header}"
"#ifndef MarlinVersion_h
#define MarlinVersion_h 1

// version macros
#define MARLIN_MAJOR_VERSION ${Marlin_MAJOR_VERSION}
#define MARLIN_MINOR_VERSION ${Marlin_MINOR_VERSION}
#define MARLIN_PATCH_LEVEL ${Marlin_PATCH_LEVEL}
#define MARLIN_VERSION_GE( MAJV , MINV , PLEV)  ( (  MARLIN_MAJOR_VERSION > MAJV ) || ( (MARLIN_MAJOR_VERSION==MAJV) && ( MARLIN_MINOR_VERSION > MINV ) ) ||   ( (MARLIN_MAJOR_VERSION==MAJV) && ( MARLIN_MINOR_VERSION == MINV ) && ( MARLIN_PATCH_LEVEL >= PLEV ) ) )

#endif
"
)

# include previously generated marlinversion.h
INCLUDE_DIRECTORIES( BEFORE "${PROJECT_BINARY_DIR}" )

# and install it
INSTALL( DIRECTORY "${PROJECT_BINARY_DIR}/marlin" DESTINATION include PATTERN "*~" EXCLUDE )

# include directories
INCLUDE_DIRECTORIES( BEFORE "${Marlin_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES( BEFORE "${Marlin_SOURCE_DIR}/streamlog/source/include" )

# install include files
INSTALL( DIRECTORY "${Marlin_SOURCE_DIR}/include/marlin" DESTINATION include
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )


####################################################################################################
### SOURCES
####################################################################################################

# get list of all source files in this directory
AUX_SOURCE_DIRECTORY( . lib_srcs )

# ProcessorLoader and main application should not go into the library
LIST( REMOVE_ITEM lib_srcs ./Marlin.cc ./ProcessorLoader.cc )

####################################################################################################
### LIBRARY
####################################################################################################

# create library

ADD_LIBRARY( lib_Marlin ${lib_srcs}  ${marlin_hdrs})


# link to libstreamlog
TARGET_LINK_LIBRARIES( lib_Marlin lib_streamlog  lib_LCIO)
TARGET_INCLUDE_DIRECTORIES(lib_Marlin INTERFACE
  "${Marlin_SOURCE_DIR}/include"
)

# change lib_target properties
SET_TARGET_PROPERTIES( lib_Marlin PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${Marlin_VERSION}
    SOVERSION ${Marlin_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME Marlin )

# install library
INSTALL( TARGETS lib_Marlin DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

####################################################################################################
### BINS
####################################################################################################

# check if we are on a *nix system
IF( NOT UNIX OR APPLE )
    IF( DEFINED LINK_STATIC_WHOLE_LIBS )
        MESSAGE( FATAL_ERROR "The LINK_STATIC_WHOLE_LIBS option is only avaliabe for *nix systems" )
    ENDIF()
    IF( NOT BUILD_SHARED_LIBS )
        MESSAGE( FATAL_ERROR "You can only build Marlin statically under *nix systems" )
    ENDIF()
ENDIF()

# create symbolic bin target for calling targets bin_XXX
ADD_CUSTOM_TARGET( bin_marlin )

IF( NOT MARLIN_NO_DLL )
    FIND_LIBRARY( DL_LIB NAMES ${CMAKE_DL_LIBS} dl )
    IF( NOT DL_LIB )
        MESSAGE( FATAL_ERROR "failed to find libdl.so" )
    ELSE()
        MESSAGE( STATUS "using dl library: ${DL_LIB}" )
    ENDIF()
    LINK_LIBRARIES( ${DL_LIB} )
ENDIF()

# whole archive libraries
IF( DEFINED LINK_STATIC_WHOLE_LIBS )
    SEPARATE_ARGUMENTS( LINK_STATIC_WHOLE_LIBS )
    SET( LINK_STATIC_WHOLE_LIBS "${LINK_STATIC_WHOLE_LIBS}" CACHE STRING
            "Link Marlin with these static libs" FORCE )
    MESSAGE( STATUS "Following static libraries are going to be linked to the binaries with the" )
    MESSAGE( STATUS "-Wl,--whole-archive linker flag:" )
    FOREACH( lib ${LINK_STATIC_WHOLE_LIBS} )
         MESSAGE( STATUS "${lib}" )
    ENDFOREACH()
ENDIF()


# ----- Marlin binary --------------------------------------------------------
ADD_EXECUTABLE( Marlin Marlin.cc ProcessorLoader.cc )
ADD_DEPENDENCIES( bin_marlin Marlin )

# link against library
IF( BUILD_SHARED_LIBS )
    TARGET_LINK_LIBRARIES( Marlin lib_Marlin )
ELSE()
    TARGET_LINK_LIBRARIES( Marlin -Wl,--whole-archive lib_Marlin -Wl,--no-whole-archive )
ENDIF()

# whole archive libraries
IF( DEFINED LINK_STATIC_WHOLE_LIBS )
    FOREACH( lib ${LINK_STATIC_WHOLE_LIBS} )
        TARGET_LINK_LIBRARIES( Marlin -Wl,--whole-archive ${lib} -Wl,--no-whole-archive )
    ENDFOREACH()
ENDIF()

INSTALL( TARGETS Marlin DESTINATION bin_marlin )
# ----------------------------------------------------------------------------



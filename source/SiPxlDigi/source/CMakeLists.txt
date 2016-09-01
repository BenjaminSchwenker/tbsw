#############################################################################
#
# CMAKE file for building SiPxlDigi
#
#
#############################################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# need long long for int64 for now
# FIXME: should establish wether we are on a 32bit or 64 bit machine...
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ../include )

# get list of all source files
AUX_SOURCE_DIRECTORY( src lib_srcs )

# create library with the list of library sources
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${GEAR_VERSION}
    SOVERSION ${GEAR_SOVERSION}
    # allow building static and shared libraries on the same build tree
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME gearxml )

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

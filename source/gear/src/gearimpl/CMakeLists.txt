##################################
# cmake file for building GEAR
# @author Jan Engels, DESY
##################################

# get list of all source files
AUX_SOURCE_DIRECTORY( src lib_srcs )
file(GLOB_RECURSE gear_hdrs  "${PROJECT_SOURCE_DIR}/src/s" *.h *.hh)

# create library with the list of library sources
ADD_LIBRARY( lib_GEAR ${lib_srcs} ${gear_hdrs} )


# change lib_target properties
SET_TARGET_PROPERTIES( lib_GEAR PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${GEAR_VERSION}
    SOVERSION ${GEAR_SOVERSION}
    # allow building static and shared libraries on the same build tree
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME gear )

# install library
INSTALL( TARGETS lib_GEAR DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )



##################################
# cmake file for building GEAR
# @author Jan Engels, DESY
##################################

# get list of all source files
AUX_SOURCE_DIRECTORY( src lib_srcs )

# create library with the list of library sources
ADD_LIBRARY( lib_GEAR ${lib_srcs} )

# cmake 2.4.6 cannot find the lib target
IF( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4.6 )
    # tell lib target to call libGEAR
    ADD_DEPENDENCIES( lib lib_GEAR )
ELSE()
    ADD_CUSTOM_TARGET( lib DEPENDS lib_GEAR )
ENDIF()

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


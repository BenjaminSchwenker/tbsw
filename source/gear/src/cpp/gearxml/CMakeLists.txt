####################################
# cmake file for building GEARXML
# @author Jan Engels, DESY
####################################

# project name
PROJECT( GEARXML )

# get list of all source files
AUX_SOURCE_DIRECTORY( src lib_srcs )

# create library with the list of library sources
ADD_LIBRARY( lib_GEARXML ${lib_srcs} )

# cmake 2.4.6 cannot find the lib target
IF( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4.6 )
    # tell lib target to call libGEARXML
    ADD_DEPENDENCIES( lib lib_GEARXML )
ELSE()
    ADD_CUSTOM_TARGET( lib DEPENDS lib_GEARXML )
ENDIF()

# link gear library with gear
TARGET_LINK_LIBRARIES( lib_GEARXML lib_GEAR )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_GEARXML PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${GEAR_VERSION}
    SOVERSION ${GEAR_SOVERSION}
    # allow building static and shared libraries on the same build tree
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME gearxml )

# install library
INSTALL( TARGETS lib_GEARXML DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )


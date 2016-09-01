#################################
# cmake file for building LCCD
# @author Jan Engels, DESY
#################################

####################################################################################################
### INCLUDES
####################################################################################################

# include directories
INCLUDE_DIRECTORIES( BEFORE "${PROJECT_SOURCE_DIR}/source/include" )

# install include files
INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/source/include" DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )


####################################################################################################
### SOURCES
####################################################################################################

# get list of all source files
AUX_SOURCE_DIRECTORY( src lib_srcs )

SET_SOURCE_FILES_PROPERTIES(
  ${lib_srcs} PROPERTIES
  COMPILE_FLAGS "-Wall -ansi -pedantic"
)

####################################################################################################
### LIBRARY
####################################################################################################

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME lccd )

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )


####################################################################################################
### BINS
####################################################################################################

# LCCD tests
SET( lccd_tests
    createcalfile
    readcalfile
    readdatafile
    testdbfile
)

# CondDB tests
SET( conddb_tests
    calfilefromdb
    calfiletodb
    createdbfile
    tagdbfolder
    multiblefolders
)

IF( BUILD_LCCD_TESTS )
    LIST( APPEND test_bins ${lccd_tests} )
ENDIF()
IF( BUILD_CONDDB_TESTS )
    LIST( APPEND test_bins ${conddb_tests} )
ENDIF()


IF( test_bins )
    # create symbolic bin target for calling targets bin_XXX
    ADD_CUSTOM_TARGET( bin )

    # loop through list of binaries
    FOREACH( bin ${test_bins} )
        SET_SOURCE_FILES_PROPERTIES( test/${bin}.cc PROPERTIES
            # FIXME -pedantic crashes with gcc 3.4.6
          COMPILE_FLAGS "-Wall -ansi"
        )
        # create target for binary
        ADD_EXECUTABLE( bin_${bin} test/${bin}.cc test/CalibrationConstant.cc )
        # tell bin target to call this bin_XXX target
        ADD_DEPENDENCIES( bin bin_${bin} )
        # avoid conflicts between library and binary target names
        SET_TARGET_PROPERTIES( bin_${bin} PROPERTIES OUTPUT_NAME ${bin} )
        # link against project library and DLL lib
        TARGET_LINK_LIBRARIES( bin_${bin} lib_${PROJECT_NAME} )
        # install binary commando
        INSTALL( TARGETS bin_${bin} DESTINATION test )
    ENDFOREACH()
ENDIF()


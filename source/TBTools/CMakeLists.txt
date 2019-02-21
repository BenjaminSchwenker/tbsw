#############################################################################
#
# CMAKE file for building TBTools
# 
# @author Benjamin Schwenker, Uni Göttingen
#
#############################################################################


SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

### PROJECT #################################################################

# Project name
PROJECT( TBTools )

# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


### create root dictionary #################################################

# create a new dictionary for class AlignEvent
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/include/AlignEventDict.h" )
    MESSAGE( STATUS "Creating root dictionary for class AlignEvent" )
    
    EXEC_PROGRAM( ${ROOT_rootcint_CMD} "${PROJECT_SOURCE_DIR}/include"
            ARGS -f AlignEventDict.C 
            ARGS -c AlignEvent.h LinkDef.h
            OUTPUT_VARIABLE out_tmp
            RETURN_VALUE out_ret )
    IF( out_ret )
        MESSAGE( FATAL_ERROR "Error executing rootcint!!" )
    ENDIF()
    
ENDIF()

# move dictionary file to other sources 
IF( EXISTS "${PROJECT_SOURCE_DIR}/include/AlignEventDict.C" )
    MESSAGE( STATUS "Created dictionary for class AlignEvent" )
    
    EXEC_PROGRAM( mv "${PROJECT_SOURCE_DIR}/include"
            ARGS AlignEventDict.C 
            ARGS ${PROJECT_SOURCE_DIR}/src
            OUTPUT_VARIABLE out_tmp
            RETURN_VALUE out_ret )
    IF( out_ret )
        MESSAGE( FATAL_ERROR "Error moving dictionary!!" )
    ENDIF()
    
ENDIF()

### INCLUDE, SOURCE, LIB ####################################################

# Include directories

INCLUDE_DIRECTORIES( include ${Eigen_INCLUDE_PATH})
INCLUDE_DIRECTORIES(SYSTEM ${ROOT_INCLUDE_DIRS})
# Set list of sources
AUX_SOURCE_DIRECTORY( src lib_srcs)
file(GLOB_RECURSE lib_hdrs RELATIVE "${PROJECT_SOURCE_DIR}" *.h *.hh)

# Create library with the list of library sources
ADD_LIBRARY( ${PROJECT_NAME} ${lib_srcs} ${lib_hdrs} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} lib_LCIO lib_Marlin ${ROOT_LIBRARIES})

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} INTERFACE
  "${PROJECT_SOURCE_DIR}/include"
)

# Install libraries
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${PROJECT_SOURCE_DIR}/include/AlignEventDict_rdict.pcm
                ${LIBRARY_OUTPUT_PATH}/AlignEventDict_rdict.pcm)

### FINAL STATUS ############################################################

# Display status message for important variables
if (VERBOSE)
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------- ${PROJECT_NAME} --------------------------------------" )
MESSAGE( STATUS "${PROJECT_NAME}_VERSION = ${${PROJECT_NAME}_VERSION}" )
ENDIF()

#############################################################################
#
# CMAKE file for building X0Tools
# 
# @author Benjamin Schwenker, Uni Göttingen
#
#############################################################################


SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

### PROJECT #################################################################

# Project name
PROJECT( X0Tools )

# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )



# Project dependencies
IF( PROJECT_DEPENDS )
    SET( ${PROJECT_NAME}_DEPENDS "${PROJECT_DEPENDS}" )
ELSE()
    SET( ${PROJECT_NAME}_DEPENDS "ROOT")
ENDIF()    

### ADVANCED SETTINGS #######################################################

# Library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )


####################################################################################################
### BINS
####################################################################################################


# Create symbolic target for calling targets bin_XXX
ADD_CUSTOM_TARGET( bin_x0 )
 
# ==== a macro for adding binaries =====
MACRO( ADD_X0TOOLS_BIN file )
  ADD_EXECUTABLE( bin_${file} ./source/src/${file}.C )
  ADD_DEPENDENCIES( bin_x0 bin_${file} )
  SET_TARGET_PROPERTIES( bin_${file} PROPERTIES OUTPUT_NAME ${file} )
  TARGET_LINK_LIBRARIES( bin_${file} ${ROOT_LIBRARIES} )
  INSTALL( TARGETS bin_${file} DESTINATION bin_x0 )
ENDMACRO()
#===================================================================

ADD_X0TOOLS_BIN( calibrationfit ) 
ADD_X0TOOLS_BIN( DrawBoxes ) 
ADD_X0TOOLS_BIN( MergeImages ) 
ADD_X0TOOLS_BIN( x0imaging ) 

### FINAL STATUS ############################################################

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "------------------------- ${PROJECT_NAME} --------------------------------------" )
MESSAGE( STATUS "${PROJECT_NAME}_VERSION = ${${PROJECT_NAME}_VERSION}" )






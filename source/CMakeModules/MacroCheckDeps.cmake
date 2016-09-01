############################################################
# cmake macro for checking dependencies
#
# Macro for looping through variables:
#   ${PROJECT_NAME}_DEPENDS
#   BUILD_WITH
#   LINK_WITH
# and load the project required/optional dependencies
# by calling the LOAD_PACKAGE macro for each of them with
# the appropriate arguments
#
# see MacroLoadPackage.cmake for more details.
#
# @author Jan Engels, DESY
############################################################

MACRO( CHECK_DEPS )

    # load macro
    MESSAGE( STATUS "-------------------------------------------------------------------------------" )
    MESSAGE( STATUS "Change a module with: cmake -D<ModuleName>_HOME=<Path_to_Module>" )
    MESSAGE( STATUS )

    # load macro
    IF( NOT DEFINED CMAKE_MODULE_PATH )
        MESSAGE( STATUS "Warning: CMAKE_MODULE_PATH not set! You can set it with: "
            "-DCMAKE_MODULE_PATH=\"/path_to_ilcsoft/CMakeModules;/foo/bar\"" )
    ENDIF()
    INCLUDE( "MacroLoadPackage" )

    # project dependencies
    IF( DEFINED ${PROJECT_NAME}_DEPENDS )
        SEPARATE_ARGUMENTS( ${PROJECT_NAME}_DEPENDS )
        MARK_AS_ADVANCED( ${PROJECT_NAME}_DEPENDS )
        FOREACH( req_pkg ${${PROJECT_NAME}_DEPENDS} )
            LOAD_PACKAGE( ${req_pkg} REQUIRED )
        ENDFOREACH()
        SET( ${PROJECT_NAME}_DEPENDS "${${PROJECT_NAME}_DEPENDS}" CACHE STRING
            "${PROJECT_NAME} dependencies" FORCE )
    ENDIF()

    # user defined dependencies
    IF( DEFINED BUILD_WITH )
        SEPARATE_ARGUMENTS( BUILD_WITH )
        MARK_AS_ADVANCED( BUILD_WITH )
        FOREACH( opt_pkg ${BUILD_WITH} )
            LOAD_PACKAGE( ${opt_pkg} REQUIRED )
        ENDFOREACH()
        SET( BUILD_WITH "${BUILD_WITH}" CACHE STRING
            "Build ${PROJECT_NAME} with these optional packages" FORCE )
    ENDIF()

    # user defined dependencies
    IF( DEFINED LINK_WITH )
        SEPARATE_ARGUMENTS( LINK_WITH )
        MARK_AS_ADVANCED( LINK_WITH )
        FOREACH( lnk_pkg ${LINK_WITH} )
            LOAD_PACKAGE( ${lnk_pkg} REQUIRED LINK_ONLY )
        ENDFOREACH()
        SET( LINK_WITH "${LINK_WITH}" CACHE STRING
            "Link ${PROJECT_NAME} with these optional packages" FORCE )
    ENDIF()

    MESSAGE( STATUS "-------------------------------------------------------------------------------" )
ENDMACRO( CHECK_DEPS )

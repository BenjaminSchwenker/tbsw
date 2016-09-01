#######################################################################
# cmake macro for loading ILC Software Packages
#
# LOAD_PACKAGE( Package_Name [QUIET] [REQUIRED] [LINK_ONLY] [NOACTION] )
#
# QUIET     no error messages unless package is required and not found
# REQUIRED  throws a fatal error if package is not found
# LINK_ONLY just links against the package libraries
# NOACTION  use this if you don't want the macro to automatically
#           incude dirs, add definitions and link against libraries.
#           If you use this option the variables set by the macro,
#           or better, set by the packages are:
#           <Package_Name>_INCLUDE_DIRS
#           <Package_Name>_LIBRARIES
#           <Package_Name>_DEFINITIONS
#
#
# WARNING:  The macro assumes that a <Package_Name>_HOME variable is
#           set and that it points to the location of the package
#           you want to load!! If this variable is not set the macro
#           will send an error or abort the configuration depending
#           on the REQUIRED argument.
#
# @author Jan Engels, DESY
#
#######################################################################

MACRO( LOAD_PACKAGE Name )

    IF( NOT ${Name}_CHECKED )
        SET( ${Name}_CHECKED TRUE INTERNAL "Resolve cyclic dependencies" )

        IF( NOT ${Name} STREQUAL ${PROJECT_NAME} )
            # initialize variables
            SET( args "" )
            SET( quiet FALSE )
            SET( required FALSE )
            SET( link_only FALSE )
            SET( no_action FALSE )

            # process arguments
            FOREACH( arg ${ARGN} )
                IF( ${arg} MATCHES "QUIET" )
                    SET( quiet TRUE )
                    SET( args ${args} ${arg} )
                ENDIF()
                IF( ${arg} MATCHES "REQUIRED" )
                    SET( required TRUE )
                    SET( args ${args} ${arg} )
                ENDIF()
                IF( ${arg} MATCHES "LINK_ONLY" )
                    SET( link_only TRUE )
                ENDIF()
                IF( ${arg} MATCHES "NOACTION" )
                    SET( no_action TRUE )
                ENDIF()
            ENDFOREACH()

            IF( NOT DEFINED ${Name}_HOME )
                IF( NOT quiet )
                    MESSAGE( STATUS "Check for ${Name}:" )
                    MESSAGE( STATUS "Check for ${Name}: ${Name}_HOME not set!! "
                        "Use: cmake -D${Name}_HOME=<path_to_${Name}>" )
                ENDIF()
                FIND_PACKAGE( ${Name} ${args} )
                IF( ${Name}_FOUND )
                    MESSAGE( STATUS "Check for ${Name}: ${Name} autodetected successfully!" )
                    IF( NOT quiet )
                        MESSAGE( STATUS "${Name}_INCLUDE_DIRS: ${${Name}_INCLUDE_DIRS}" )
                        MESSAGE( STATUS "${Name}_LIBRARIES: ${${Name}_LIBRARIES}" )
                        MESSAGE( STATUS "${Name}_DEFINITIONS: ${${Name}_DEFINITIONS}" )
                    ENDIF()
                ELSE()
                    IF( required )
                        MESSAGE( FATAL_ERROR "Check for ${Name}: failed to autodetect ${Name}!!" )
                    ELSE()
                        MESSAGE( STATUS "Check for ${Name}: failed to autodetect ${Name}!!" )
                    ENDIF()
                ENDIF()
            ELSE()
                # create a cache entry for the Package_HOME variable
                # this is needed because command line defined variables
                # (-DVariable=Value) are not cached by default
                SET( ${Name}_HOME ${${Name}_HOME} CACHE PATH "Path to ${Name}" FORCE )


                # if a Config.cmake exists use it!
                IF( EXISTS ${${Name}_HOME}/${Name}Config.cmake )
                    
                    MESSAGE( STATUS "Check for ${Name}: ${${Name}_HOME} -- ${Name}Config.cmake found..." )
                    
                    # this variable needs to be set for cmake to read the config file
                    SET( ${Name}_DIR ${${Name}_HOME} CACHE PATH "Path to ${Name}" FORCE )
                    MARK_AS_ADVANCED( ${Name}_DIR )
                    
                    # as opposed to cmake Find modules the
                    # Config.cmake files are not aware of this variable
                    SET( ${Name}_FIND_REQUIRED ${required} )
                    
                    FIND_PACKAGE( ${Name} ${args} NO_MODULE )
                ELSE()
                    MESSAGE( STATUS "Check for ${Name}: ${${Name}_HOME}" )
                    FIND_PACKAGE( ${Name} ${args} )
                ENDIF()
            ENDIF()
            # include directories
            IF( ${Name}_INCLUDE_DIRS )
                IF( NOT link_only AND NOT no_action )
                    INCLUDE_DIRECTORIES( ${${Name}_INCLUDE_DIRS} )
                ENDIF()
            ELSE()
                IF( ${Name}_INCLUDE_DIR )
                    # set <PKG>_INCLUDE_DIRS if not set
                    SET( ${Name}_INCLUDE_DIRS ${${Name}_INCLUDE_DIR} )
                    IF( NOT link_only AND NOT no_action )
                        INCLUDE_DIRECTORIES( ${${Name}_INCLUDE_DIR} )
                    ENDIF()
                ENDIF()
            ENDIF()

            # set library dependencies
            IF( ${Name}_LIBRARIES )
                SEPARATE_ARGUMENTS( ${Name}_LIBRARIES )
                LIST( APPEND ${PROJECT_NAME}_LIBRARIES ${${Name}_LIBRARIES} )
                IF( NOT no_action )
                    LINK_LIBRARIES( ${${Name}_LIBRARIES} )
                ENDIF()
            ENDIF()

            # add definitions
            IF( ${Name}_DEFINITIONS )
                #FIXME not viable to do this, yet!
                #SEPARATE_ARGUMENTS( ${Name}_DEFINITIONS )
                ##remove duplicates
                #FOREACH( new_def ${${Name}_DEFINITIONS})
                #    SET( def_found FALSE )
                #    FOREACH( def ${${PROJECT_NAME}_DEFINITIONS})
                #        IF( "${def}" STREQUAL "${new_def}" )
                #            SET( def_found TRUE )
                #        ENDIF()
                #    ENDFOREACH()
                #    IF( NOT def_found )
                #        LIST( APPEND ${PROJECT_NAME}_DEFINITIONS ${new_def} )
                #    ENDIF()
                #ENDFOREACH()
                    
                #LIST( APPEND ${PROJECT_NAME}_DEFINITIONS ${${Name}_DEFINITIONS} )
                IF( NOT no_action )
                    ADD_DEFINITIONS( ${${Name}_DEFINITIONS} )
                ENDIF()
            ENDIF()
        ENDIF()
    ENDIF()

ENDMACRO( LOAD_PACKAGE )

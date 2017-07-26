#################################
# cmake file for building LCIO
# @author Frank Gaede, DESY
# @author Jan Engels, DESY
#################################

INCLUDE_DIRECTORIES (
    "${LCIO_SOURCE_DIR}/sio/include"
    "${LCIO_SOURCE_DIR}/src/cpp/include"
)



SET( LCIO_IMPL_SRCS
  ./src/IMPL/AccessChecked.cc
  ./src/IMPL/CalorimeterHitImpl.cc
  ./src/IMPL/ClusterImpl.cc
  ./src/IMPL/LCCollectionVec.cc
  ./src/IMPL/LCEventImpl.cc
  ./src/IMPL/LCFlagImpl.cc
  ./src/IMPL/LCGenericObjectImpl.cc
  ./src/IMPL/LCIO.cc
  ./src/IMPL/LCIOExceptionHandler.cc
  ./src/IMPL/LCParametersImpl.cc
  ./src/IMPL/LCRunHeaderImpl.cc
  ./src/IMPL/MCParticleImpl.cc
  ./src/IMPL/ParticleIDImpl.cc
  ./src/IMPL/RawCalorimeterHitImpl.cc
  ./src/IMPL/ReconstructedParticleImpl.cc
  ./src/IMPL/SimCalorimeterHitImpl.cc
  ./src/IMPL/SimTrackerHitImpl.cc
  ./src/IMPL/TPCHitImpl.cc
  ./src/IMPL/TrackImpl.cc
  ./src/IMPL/TrackerDataImpl.cc
  ./src/IMPL/TrackerHitImpl.cc
  ./src/IMPL/TrackerPulseImpl.cc
  ./src/IMPL/TrackerRawDataImpl.cc
  ./src/IMPL/VertexImpl.cc
)

SET( LCIO_IOIMPL_SRCS
  ./src/IOIMPL/LCFactory.cc
)

SET( LCIO_SIO_SRCS
  ./src/SIO/LCSIO.cc
  ./src/SIO/SIOCalHitHandler.cc
  ./src/SIO/SIOClusterHandler.cc
  ./src/SIO/SIOCollectionHandler.cc
  ./src/SIO/SIOEventHandler.cc
  ./src/SIO/SIOFloatVecHandler.cc
  ./src/SIO/SIOHandlerMgr.cc
  ./src/SIO/SIOIntVecHandler.cc
  ./src/SIO/SIOLCGenericObjectHandler.cc
  ./src/SIO/SIOLCParameters.cc
  ./src/SIO/SIOLCRelationHandler.cc
  ./src/SIO/SIOObjectHandler.cc
  ./src/SIO/SIOParticleHandler.cc
  ./src/SIO/SIORawCalHitHandler.cc
  ./src/SIO/SIOReconstructedParticleHandler.cc
  ./src/SIO/SIORunHeaderHandler.cc
  ./src/SIO/SIOSimCalHitHandler.cc
  ./src/SIO/SIOSimTrackHitHandler.cc
  ./src/SIO/SIOStrVecHandler.cc
  ./src/SIO/SIOTPCHitHandler.cc
  ./src/SIO/SIOTrackHandler.cc
  ./src/SIO/SIOTrackerDataHandler.cc
  ./src/SIO/SIOTrackerHitHandler.cc
  ./src/SIO/SIOTrackerPulseHandler.cc
  ./src/SIO/SIOTrackerRawDataHandler.cc
  ./src/SIO/SIOUnpack.cc
  ./src/SIO/SIOVertexHandler.cc
  ./src/SIO/SIOWriter.cc
  ./src/SIO/SIOReader.cc
)

SET( LCIO_UTIL_SRCS
  ./src/UTIL/BitField64.cc
  ./src/UTIL/CellIDDecoder.cc
  ./src/UTIL/IndexMap.cc
  ./src/UTIL/LCRelationNavigator.cc
  ./src/UTIL/LCSplitWriter.cc
  ./src/UTIL/LCStdHepRdr.cc
  ./src/UTIL/LCTime.cc
  ./src/UTIL/Operators.cc
  ./src/UTIL/lStdHep.cc
  ./src/UTIL/lXDR.cc
  ./src/UTIL/LCTOOLS.cc
  ./src/UTIL/CellIDEncoder.cc
  ./src/UTIL/CollectionParameterMap.cc
  ./src/UTIL/PIDHandler.cc
)

#------ require proper C++ -------
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )


ADD_CUSTOM_TARGET( lib DEPENDS lib_LCIO lib_SIO )

ADD_LIBRARY( lib_LCIO 
  ${LCIO_IMPL_SRCS}
  ${LCIO_IOIMPL_SRCS}
  ${LCIO_UTIL_SRCS}
  ${LCIO_SIO_SRCS}
)

SET_TARGET_PROPERTIES( lib_LCIO PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${LCIO_VERSION}
    SOVERSION ${LCIO_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME lcio
)

TARGET_LINK_LIBRARIES( lib_LCIO lib_SIO )


# create symbolic target for calling targets bin_XXX
ADD_CUSTOM_TARGET( bin )

# ==== a macro for adding binaries from the EXAMPLE directory =====
MACRO( ADD_LCIO_BIN file )
  ADD_EXECUTABLE( bin_${file} ./src/EXAMPLE/${file}.cc )
  ADD_DEPENDENCIES( bin bin_${file} )
#  MESSAGE( STATUS " in ADD_LCIO_BIN : " ${file} " " ./src/EXAMPLE/${file}.cc )
  SET_TARGET_PROPERTIES( bin_${file} PROPERTIES OUTPUT_NAME ${file} )
  TARGET_LINK_LIBRARIES( bin_${file} lib_LCIO )
  INSTALL( TARGETS bin_${file} DESTINATION bin )
ENDMACRO()
#===================================================================

ADD_LCIO_BIN( dumpevent ) 
ADD_LCIO_BIN( stdhepjob ) 

#------ create a root dictionary -----------------------------
IF( BUILD_ROOTDICT )
#-----------------------------------------------------------------

  INCLUDE( ${PROJECT_SOURCE_DIR}/FindROOT.cmake )


  IF( ROOT_FOUND )
    
    
    SET( ENV{ROOTSYS} "${ROOT_HOME}" )
    SET( ENV{PATH} $ENV{ROOTSYS}/bin:$ENV{PATH} )
    SET( ENV{LD_LIBRARY_PATH} $ENV{ROOTSYS}/lib:$ENV{LD_LIBRARY_PATH} )
    
    
    #TODO test ROOT_FOUND and ${ROOT_HOME}/bin/rootcint
    
    # list of headers
    FILE( GLOB_RECURSE lcio_headers include/*.h )
    
    # list of headers to exclude
    FILE( GLOB_RECURSE exclude_headers include/empty_ignore.h include/SIO/*.h include/CPPFORT/*.h include/UTIL/CellIDDecoder.h)
    
    # remove items to exclude
    LIST( REMOVE_ITEM lcio_headers ${exclude_headers} )
    
    SET( DICT_OUTPUT_DIR "${PROJECT_BINARY_DIR}/dict" )
    
    # create output directory
#    ADD_CUSTOM_COMMAND(
#      OUTPUT  "${DICT_OUTPUT_DIR}"
#      COMMAND "${CMAKE_COMMAND}"
#      ARGS -E make_directory "${DICT_OUTPUT_DIR}"
#      COMMENT "create directory ${DICT_OUTPUT_DIR}"
#      )
    file(MAKE_DIRECTORY ${DICT_OUTPUT_DIR})



    SET( dict_lib_srcs )
    
    # FIXME variables must be set by FindROOT.cmake
    #SET( cint_script "${PROJECT_BINARY_DIR}/rootcint.cmake" )
    #FILE( WRITE "${cint_script}" "# little script to run rootcint\n" )
    #FILE( APPEND "${cint_script}" "SET( ENV{ROOTSYS} ${ROOT_HOME} )\n" )
    #FILE( APPEND "${cint_script}" "SET( ENV{PATH} \$ENV{ROOTSYS}/bin:\$ENV{PATH} )\n" )
    #FILE( APPEND "${cint_script}" "SET( ENV{LD_LIBRARY_PATH} \$ENV{ROOTSYS}/lib:\$ENV{LD_LIBRARY_PATH} )\n" )
    #FILE( APPEND "${cint_script}" "EXEC_PROGRAM( rootcint -f \${dict_src_file} -c -DHANDLE_LCIO_EXCEPTIONS=IGNORED_FOR_CINT -Iinclude \${lcio_hdr_rel_path} \${dict_hdr_file} )\n" )
    
    
    FOREACH( lcio_hdr ${lcio_headers} )
      
      # get the relative path
      #FILE( RELATIVE_PATH lcio_hdr_rel_path "${PROJECT_SOURCE_DIR}" ${lcio_hdr} )
      #STRING( REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" lcio_hdr_rel_path ${lcio_hdr} )
      FILE( RELATIVE_PATH lcio_hdr_rel_path "${CMAKE_CURRENT_SOURCE_DIR}" ${lcio_hdr} )
      
      # generate a basename from the header + relative path, e.g.:
      # include/EVENT/TrackerPulse.h  -->>  include_EVENT_TrackerPulse
      STRING( REPLACE "/" "_" dict_basename ${lcio_hdr_rel_path} )
      STRING( REPLACE ".h" "" dict_basename ${dict_basename} )
      
      
      SET( dict_hdr "${dict_basename}_Linkdef.h" )
      SET( dict_hdr_file "${DICT_OUTPUT_DIR}/${dict_hdr}" )
      
      SET( dict_src "${dict_basename}_dict.cxx" )
      SET( dict_src_file "${DICT_OUTPUT_DIR}/${dict_src}" )
      
      # generate dict header file
      ADD_CUSTOM_COMMAND(
	OUTPUT  ${dict_hdr_file}
	COMMAND ${CMAKE_COMMAND}
	ARGS -E echo "\\#pragma link C++ defined_in \\\"${lcio_hdr_rel_path}\\\"\; >${dict_hdr}"
	#ARGS -E echo "\\#pragma link C++ defined_in ${lcio_hdr_rel_path}\; >${dict_hdr}"
	WORKING_DIRECTORY "${DICT_OUTPUT_DIR}"
	DEPENDS ${lcio_hdr}  # ${DICT_OUTPUT_DIR}
	COMMENT "generating ${dict_hdr}"
	)
      
      # generate dict source file
      ADD_CUSTOM_COMMAND(
	OUTPUT  ${dict_src_file}
	#COMMAND ${ROOT_HOME}/bin/rootcint
	#COMMAND $ENV{ROOTSYS}/bin/rootcint
	COMMAND rootcint
	ARGS -f ${dict_src_file} -c -DHANDLE_LCIO_EXCEPTIONS=IGNORED_FOR_CINT -Iinclude ${lcio_hdr_rel_path} ${dict_hdr_file}
	WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
	DEPENDS ${dict_hdr_file}
	COMMENT "generating ${dict_src}"
	)
      
      LIST( APPEND dict_lib_srcs ${dict_src_file} )
    ENDFOREACH()
    
    INCLUDE_DIRECTORIES( "${DICT_OUTPUT_DIR}" )
    INCLUDE_DIRECTORIES( "${ROOT_INCLUDE_DIRS}" )
    # FIXME next line should not be needed
    # xxx_dict.h files should have:
    # #include "IOIMPL/TrackerRawDataIOImpl.h"
    # instead of:
    # #include "include/IOIMPL/TrackerRawDataIOImpl.h"
    INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}" )
    
    ADD_LIBRARY( lib_LCIODICT ${dict_lib_srcs} ${LCIO_RIO_SRCS} )
    INSTALL( TARGETS lib_LCIODICT DESTINATION lib )
 

   IF( APPLE )
     SET( ROOT_DICT_LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load")
#     ADD_LIBRARY( lib_LCIO ${LCIO_RIO_SRCS} )
   ELSE( APPLE )
     SET( ROOT_DICT_LINK_FLAGS " " )
   ENDIF( APPLE )
    
    SET_TARGET_PROPERTIES( lib_LCIODICT PROPERTIES
      # create *nix style library versions + symbolic links
      VERSION ${LCIO_VERSION}
      SOVERSION ${LCIO_SOVERSION}
      # allow creating static and shared libs without conflicts
      CLEAN_DIRECT_OUTPUT 1
      # avoid conflicts between library and binary target names
      OUTPUT_NAME lcioDict
      LINKER_FLAGS ${ROOT_DICT_LINK_FLAGS}
      )
    
    ADD_CUSTOM_TARGET( dict DEPENDS  lib_LCIODICT )

    TARGET_LINK_LIBRARIES( lib_LCIODICT ${ROOT_LIBRARIES} lib_LCIO)


  ELSE( ROOT_FOUND ) 

     MESSAGE( STATUS " NO ROOT found :-( " ) 

     
  ENDIF( ROOT_FOUND )
ENDIF( BUILD_ROOTDICT )


#------------- install the library and the (public) header files -------

INSTALL( TARGETS lib_LCIO DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )



INSTALL(DIRECTORY ${LCIO_SOURCE_DIR}/src/cpp/include/ DESTINATION include 
  PATTERN "*~" EXCLUDE
  #PATTERN "*empty_ignore.h" EXCLUDE
  PATTERN "*SIO*" EXCLUDE
  PATTERN "*CPPFORT*" EXCLUDE
  PATTERN "*CVS*" EXCLUDE)



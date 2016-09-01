#############################################
# cmake file for building GEAR test binaries
# @author Jan Engels, DESY
#############################################

# mergeXML
ADD_EXECUTABLE( mergeXML mergeXML.cc )
TARGET_LINK_LIBRARIES( mergeXML lib_GEAR lib_GEARXML )
INSTALL( TARGETS mergeXML DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

# printParameters
ADD_EXECUTABLE( printParameters printParameters.cc )
TARGET_LINK_LIBRARIES( printParameters lib_GEAR lib_GEARXML )
INSTALL( TARGETS printParameters DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

# test binaries
SET( test_bins
    testgear
    testtpcproto
    testVTXgear
    testSiPlanesgear
)

ADD_CUSTOM_TARGET( tests )

FOREACH( bin ${test_bins} )
    ADD_EXECUTABLE( ${bin} EXCLUDE_FROM_ALL ${bin}.cc )
    ADD_DEPENDENCIES( tests ${bin} )
    TARGET_LINK_LIBRARIES( ${bin} lib_GEAR lib_GEARXML )
ENDFOREACH()


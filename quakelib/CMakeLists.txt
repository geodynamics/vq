# Create a library called "quakelib" which includes the appropriate source files.
# The extension is already found. Any number of sources could be listed here. 

SET(QUAKELIB_SOURCES
    ${QUAKELIB_SOURCE_DIR}/QuakeLib.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibUtil.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibOkada.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibCondition.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibSectionParams.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibElement.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibEvent.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibFriction.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibGeometry.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibMetadata.cpp
	${QUAKELIB_SOURCE_DIR}/QuakeLibOkada.cpp
	${QUAKELIB_SOURCE_DIR}/QuakeLibUtil.cpp)

# Build the QuakeLib library
add_library (quakelib ${QUAKELIB_SOURCES})

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(quakelib PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

SET(QUAKELIB_HEADERS
    ${QUAKELIB_SOURCE_DIR}/QuakeLib.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibUtil.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibOkada.h)

INSTALL(FILES ${QUAKELIB_SOURCES} DESTINATION include/quakelib/)

INSTALL(TARGETS quakelib LIBRARY DESTINATION lib/quakelib/ COMPONENT libraries)

# Go into the python subdirectory
add_subdirectory(python)

# If we have Python, run a set of tests with Python based scripts
# TODO: change this to Python check
#if (DEFINED SWIG_FOUND)
#    add_test(NAME CondUnitTest COMMAND ${QUAKELIB_TEST_DIR}/CondUnitTest.py)
#    set_property(TEST CondUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME FricUnitTest COMMAND ${QUAKELIB_TEST_DIR}/FricUnitTest.py)
#    set_property(TEST FricUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME GreenUnitTest COMMAND ${QUAKELIB_TEST_DIR}/GreenUnitTest.py)
#    set_property(TEST GreenUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME OctreeTest COMMAND ${QUAKELIB_TEST_DIR}/OctreeTest.py)
#    set_property(TEST OctreeTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME UtilUnitTest COMMAND ${QUAKELIB_TEST_DIR}/UtilUnitTest.py)
#    set_property(TEST UtilUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME EventUnitTest COMMAND ${QUAKELIB_TEST_DIR}/EventUnitTest.py)
#    set_property(TEST EventUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME GeomUnitTest COMMAND ${QUAKELIB_TEST_DIR}/GeomUnitTest.py)
#    set_property(TEST GeomUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME MetadataUnitTest COMMAND ${QUAKELIB_TEST_DIR}/MetadataUnitTest.py)
#    set_property(TEST MetadataUnitTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#    add_test(NAME RectBoundTest COMMAND ${QUAKELIB_TEST_DIR}/RectBoundTest.py)
#    set_property(TEST RectBoundTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/quakelib/python/")
#endif(DEFINED SWIG_FOUND)


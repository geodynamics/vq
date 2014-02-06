# Create a library called "quakelib" which includes the appropriate source files.

PROJECT (QuakeLib)
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Check for common include files
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("float.h" QUAKELIB_HAVE_FLOAT_H)
CHECK_INCLUDE_FILES ("limits.h" QUAKELIB_HAVE_LIMITS_H)
CHECK_INCLUDE_FILES ("math.h" QUAKELIB_HAVE_MATH_H)
CHECK_INCLUDE_FILES ("stdlib.h" QUAKELIB_HAVE_STDLIB_H)
CHECK_INCLUDE_FILES ("string.h" QUAKELIB_HAVE_STRING_H)

# Check if Doxygen is installed
FIND_PACKAGE(Doxygen)

# Check for SWIG and Python installation
# TODO: add IF statements for INCLUDE paths
FIND_PACKAGE(SWIG)
IF(SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
ENDIF(SWIG_FOUND)

FIND_PACKAGE(PythonLibs)
IF(DEFINED PYTHONLIBS_FOUND)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF(DEFINED PYTHONLIBS_FOUND)

# Set up directory definitions
SET(QUAKELIB_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/)
SET(QUAKELIB_PYTHON_DIR ${CMAKE_SOURCE_DIR}/python/)
SET(QUAKELIB_TEST_DIR ${CMAKE_SOURCE_DIR}/test/)
INCLUDE_DIRECTORIES(${QUAKELIB_SOURCE_DIR})

# Generate the config.h file and add the include path for it
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/quakelib_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/quakelib_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(QUAKELIB_SOURCES
    ${QUAKELIB_SOURCE_DIR}/QuakeLibElement.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibEQSim.cpp
    ${QUAKELIB_SOURCE_DIR}/QuakeLibIO.cpp
	${QUAKELIB_SOURCE_DIR}/QuakeLibOkada.cpp
	${QUAKELIB_SOURCE_DIR}/QuakeLibUtil.cpp
    )

# Build the QuakeLib library
add_library (quakelib ${QUAKELIB_SOURCES})

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(quakelib PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

SET(QUAKELIB_HEADERS
    ${QUAKELIB_SOURCE_DIR}/QuakeLib.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibEQSim.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibIO.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibOkada.h
    ${QUAKELIB_SOURCE_DIR}/QuakeLibUtil.h
    )

INSTALL(FILES ${QUAKELIB_SOURCES} DESTINATION include/quakelib/)

INSTALL(TARGETS quakelib
     RUNTIME DESTINATION bin/quakelib/ COMPONENT libraries
     LIBRARY DESTINATION lib/quakelib/ COMPONENT libraries
     ARCHIVE DESTINATION lib/static/quakelib/ COMPONENT libraries)

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


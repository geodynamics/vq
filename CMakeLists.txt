# The name of our project is "VC". CMakeLists files in this project can 
# refer to the root source directory of the project as ${VC_SOURCE_DIR} and 
# to the root binary directory of the project as ${VC_BINARY_DIR}. 
cmake_minimum_required (VERSION 2.6) 
project (VC) 

# Enable testing
#INCLUDE (CTest)
# Check for common include files
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("sys/time.h" VC_HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES ("float.h" VC_HAVE_FLOAT_H)
CHECK_INCLUDE_FILES ("limits.h" VC_HAVE_LIMITS_H)
CHECK_INCLUDE_FILES ("math.h" VC_HAVE_MATH_H)
CHECK_INCLUDE_FILES ("stdlib.h" VC_HAVE_STDLIB_H)
CHECK_INCLUDE_FILES ("string.h" VC_HAVE_STRING_H)
CHECK_INCLUDE_FILES ("unistd.h" VC_HAVE_UNISTD_H)

SET(QUAKELIB_SOURCE_DIR ${VC_SOURCE_DIR}/quakelib/src/)
SET(QUAKELIB_PYTHON_DIR ${VC_SOURCE_DIR}/quakelib/python/)
#SET(QUAKELIB_TEST_DIR ${VC_SOURCE_DIR}/quakelib/test/)
#SET(VC_EXAMPLE_DIR ${VC_SOURCE_DIR}/examples/)

# Check for optional libraries
FIND_PACKAGE(MPI)
IF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)
    ADD_DEFINITIONS(${MPI_C_COMPILE_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
ENDIF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)

FIND_PACKAGE(OpenMP)

FIND_PACKAGE(HDF5 COMPONENTS C HL)
IF(DEFINED HDF5_FOUND)
    # If HDF5 requires parallel but MPI wasn't found, turn off HDF5
    IF(DEFINED HDF5_IS_PARALLEL AND NOT MPI_C_FOUND)
        SET(HDF5_FOUND)
    ENDIF(DEFINED HDF5_IS_PARALLEL AND NOT MPI_C_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF(DEFINED HDF5_FOUND)

# TODO: add IF statements for INCLUDE paths
FIND_PACKAGE(SWIG)
IF(DEFINED SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
ENDIF(DEFINED SWIG_FOUND)

FIND_PACKAGE(PythonLibs)
IF(DEFINED PYTHONLIBS_FOUND)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF(DEFINED PYTHONLIBS_FOUND)

# Create the config.h file and make sure everyone can find it
CONFIGURE_FILE(${VC_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Add common directories to the include search path
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${QUAKELIB_SOURCE_DIR})

# Recurse into the "quakelib" and "src" subdirectories.
ADD_SUBDIRECTORY(quakelib) 
ADD_SUBDIRECTORY(src)
#ADD_SUBDIRECTORY(examples)


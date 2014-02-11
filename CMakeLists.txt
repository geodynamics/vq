# The name of our project is "VC". CMakeLists files in this project can 
# refer to the root source directory of the project as ${VC_SOURCE_DIR} and 
# to the root binary directory of the project as ${VC_BINARY_DIR}. 
cmake_minimum_required (VERSION 2.6) 
project (VC) 

# Enable testing
INCLUDE (CTest)

# Check for common include files
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("sys/time.h" VC_HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES ("stdlib.h" VC_HAVE_STDLIB_H)
CHECK_INCLUDE_FILES ("unistd.h" VC_HAVE_UNISTD_H)

# Check for functions
INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS (sleep VC_HAVE_SLEEP_FUNC)
CHECK_FUNCTION_EXISTS (usleep VC_HAVE_USLEEP_FUNC)

SET( ENV{CMAKE_OSX_ARCHITECTURES} x86_64 )

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

# Create the config.h file and make sure everyone can find it
CONFIGURE_FILE(${VC_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Add common directories to the include search path
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/quakelib/)

# Recurse into the "quakelib" and "src" subdirectories.
ADD_SUBDIRECTORY(quakelib) 
ADD_SUBDIRECTORY(src)
#ADD_SUBDIRECTORY(examples)


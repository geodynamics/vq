############################################################
# USAGE
# cmake
#   generate Makefiles to compile/test VQ and QuakeLib
# cmake -G Xcode
#   generate XCode project to compile/test VQ and QuakeLib
# make -j <#>
#   compile VQ and QuakeLib with # processors
# make test OR ctest
#   run all possible tests to ensure correct functionality
# ctest -R <regex>
#   run tests that match <regex> expression
# make package_source
#   create source package for distribution
############################################################
# The name of our project is "VQ". CMakeLists files in this project can
# refer to the root source directory of the project as ${VQ_SOURCE_DIR} and
# to the root binary directory of the project as ${VQ_BINARY_DIR}.
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (VQ)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/quakelib/cmake/)

# Set the version number
SET(VQ_VERSION_MAJOR 3)
SET(VQ_VERSION_MINOR 1)
SET(VQ_VERSION_SUBMINOR 1)
SET(VQ_VERSION_STR "${VQ_VERSION_MAJOR}.${VQ_VERSION_MINOR}.${VQ_VERSION_SUBMINOR}")

# Enable CPack functionality
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Virtual Quake")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR ${VQ_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VQ_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VQ_VERSION_SUBMINOR})
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${CMAKE_SOURCE_DIR}/src;./src"
    "${CMAKE_SOURCE_DIR}/examples;./examples"
    "${CMAKE_SOURCE_DIR}/quakelib;./quakelib"
    "${CMAKE_SOURCE_DIR}/doc;./doc")
SET(CPACK_SOURCE_GENERATOR TGZ)
INCLUDE (CPack)

# Find a Git installation and the source Git revision if available
FIND_PACKAGE(Git)
IF(GIT_FOUND)
    INCLUDE(GetGitRevisionDescription)
    GET_GIT_HEAD_REVISION(VQ_GIT_REFSPEC VQ_GIT_SHA1)
ENDIF(GIT_FOUND)

# Enable testing
INCLUDE (CTest)

 # If the user doesn't specify, build the release version by default
IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF (NOT CMAKE_BUILD_TYPE)

# Check validity of different flags
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-O0 COMPILER_SUPPORTS_OPT0_FLAG)
CHECK_CXX_COMPILER_FLAG(-O3 COMPILER_SUPPORTS_OPT3_FLAG)
CHECK_CXX_COMPILER_FLAG(-g COMPILER_SUPPORTS_DEBUG_FLAG)

# Define debug vs release compiler flags
IF(COMPILER_SUPPORTS_OPT0_FLAG)
    SET(DEBUG_FLAGS "${DEBUG_FLAGS} -O0")
ENDIF(COMPILER_SUPPORTS_OPT0_FLAG)

IF(COMPILER_SUPPORTS_OPT3_FLAG)
    SET(RELEASE_FLAGS "${RELEASE_FLAGS} -O3")
    SET(RELWITHDEBINFO_FLAGS "${RELWITHDEBINFO_FLAGS} -O3")
ENDIF(COMPILER_SUPPORTS_OPT3_FLAG)

IF(COMPILER_SUPPORTS_DEBUG_FLAG)
    SET(DEBUG_FLAGS "${DEBUG_FLAGS} -g")
    SET(RELWITHDEBINFO_FLAGS "${RELWITHDEBINFO_FLAGS} -g")
ENDIF(COMPILER_SUPPORTS_DEBUG_FLAG)

SET(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
SET(CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${RELWITHDEBINFO_FLAGS})

# Check for common include files
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES ("sys/time.h" VQ_HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES ("math.h" VQ_HAVE_MATH_H)
CHECK_INCLUDE_FILES ("stdlib.h" VQ_HAVE_STDLIB_H)
CHECK_INCLUDE_FILES ("unistd.h" VQ_HAVE_UNISTD_H)
CHECK_INCLUDE_FILES ("string.h" VQ_HAVE_STRING_H)
CHECK_INCLUDE_FILES ("signal.h" VQ_HAVE_SIGNAL_H)

# Check for functions
INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS (sleep VQ_HAVE_SLEEP_FUNC)
CHECK_FUNCTION_EXISTS (usleep VQ_HAVE_USLEEP_FUNC)

SET( ENV{CMAKE_OSX_ARCHITECTURES} x86_64 )

# Check for optional libraries
# First, MPI for parallel runs
FIND_PACKAGE(MPI)
IF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)
    ADD_DEFINITIONS(${MPI_C_COMPILE_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
ENDIF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)

# OpenMP for SMP calculation
FIND_PACKAGE(OpenMP)

# HDF5 for large output
FIND_PACKAGE(HDF5 COMPONENTS C HL)
IF(DEFINED HDF5_FOUND)
    # If HDF5 requires parallel but MPI wasn't found, turn off HDF5
    IF(DEFINED HDF5_IS_PARALLEL AND NOT MPI_C_FOUND)
        SET(HDF5_FOUND)
    ENDIF(DEFINED HDF5_IS_PARALLEL AND NOT MPI_C_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF(DEFINED HDF5_FOUND)

# Create the config.h file and make sure everyone can find it
CONFIGURE_FILE(${VQ_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Add common directories to the include search path
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/quakelib/)

# Recurse into the "quakelib" and "src" subdirectories.
ADD_SUBDIRECTORY(quakelib)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)

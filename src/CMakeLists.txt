# Make sure the compiler can find include files from the quakelib library. 
INCLUDE_DIRECTORIES(${VQ_SOURCE_DIR}/quakelib/src)

# Set VQ source directories
SET(VQ_CORE_DIR ${VQ_SOURCE_DIR}/src/core)
SET(VQ_IO_DIR ${VQ_SOURCE_DIR}/src/io)
SET(VQ_MISC_DIR ${VQ_SOURCE_DIR}/src/misc)
SET(VQ_SIM_DIR ${VQ_SOURCE_DIR}/src/simulation)

# and include files from VQ itself
INCLUDE_DIRECTORIES (${VQ_CORE_DIR})
INCLUDE_DIRECTORIES (${VQ_IO_DIR})
INCLUDE_DIRECTORIES (${VQ_MISC_DIR})
INCLUDE_DIRECTORIES (${VQ_SIM_DIR})

# Make sure the linker can find the quakelib library once it is built. 
LINK_DIRECTORIES (${VQ_BINARY_DIR}/quakelib)

# Define sets of the source files
SET(VQ_MISC
    ${VQ_MISC_DIR}/ConfigFile.cpp
    ${VQ_MISC_DIR}/ConfigFile.h
    ${VQ_MISC_DIR}/GreensFunctions.cpp
    ${VQ_MISC_DIR}/GreensFunctions.h
    ${VQ_MISC_DIR}/MPIDebugOutputStream.cpp
    ${VQ_MISC_DIR}/MPIDebugOutputStream.h
    )

SET(VQ_SIMULATION
    ${VQ_SIM_DIR}/BASSAftershocks.cpp
    ${VQ_SIM_DIR}/BASSAftershocks.h
    ${VQ_SIM_DIR}/GracefulQuit.cpp
    ${VQ_SIM_DIR}/GracefulQuit.h
    ${VQ_SIM_DIR}/GreensInit.cpp
    ${VQ_SIM_DIR}/GreensInit.h
    ${VQ_SIM_DIR}/GreensKillInteraction.cpp
    ${VQ_SIM_DIR}/GreensKillInteraction.h
    ${VQ_SIM_DIR}/ProgressMonitor.cpp
    ${VQ_SIM_DIR}/ProgressMonitor.h
    ${VQ_SIM_DIR}/RunEvent.cpp
    ${VQ_SIM_DIR}/RunEvent.h
    ${VQ_SIM_DIR}/SanityCheck.cpp
    ${VQ_SIM_DIR}/SanityCheck.h
    ${VQ_SIM_DIR}/UpdateBlockStress.cpp
    ${VQ_SIM_DIR}/UpdateBlockStress.h
    )

SET(VQ_IO
    ${VQ_IO_DIR}/CheckpointFileOutput.cpp
    ${VQ_IO_DIR}/CheckpointFileOutput.h
    ${VQ_IO_DIR}/CheckpointFileParse.cpp
    ${VQ_IO_DIR}/CheckpointFileParse.h
    ${VQ_IO_DIR}/EventOutput.cpp
    ${VQ_IO_DIR}/EventOutput.h
    ${VQ_IO_DIR}/GreensFileOutput.cpp
    ${VQ_IO_DIR}/GreensFileOutput.h
    ${VQ_IO_DIR}/HDF5Data.cpp
    ${VQ_IO_DIR}/HDF5Data.h
    ${VQ_IO_DIR}/ReadModelFile.cpp
    ${VQ_IO_DIR}/ReadModelFile.h
    )

SET(VQ_CORE
    ${VQ_CORE_DIR}/SimFramework.cpp
    ${VQ_CORE_DIR}/SimFramework.h
    ${VQ_CORE_DIR}/SimTimer.cpp
    ${VQ_CORE_DIR}/SimTimer.h
    ${VQ_CORE_DIR}/Block.cpp
    ${VQ_CORE_DIR}/Block.h
    ${VQ_CORE_DIR}/Comm.cpp
    ${VQ_CORE_DIR}/Comm.h
    ${VQ_CORE_DIR}/CommPartition.cpp
    ${VQ_CORE_DIR}/CommPartition.h
    ${VQ_CORE_DIR}/Event.cpp
    ${VQ_CORE_DIR}/Event.h
    ${VQ_CORE_DIR}/InitBlocks.cpp
    ${VQ_CORE_DIR}/InitBlocks.h
    ${VQ_CORE_DIR}/Params.cpp
    ${VQ_CORE_DIR}/Params.h
    ${VQ_CORE_DIR}/SimData.cpp
    ${VQ_CORE_DIR}/SimData.h
    ${VQ_CORE_DIR}/SimDataBlocks.cpp
    ${VQ_CORE_DIR}/SimDataBlocks.h
    ${VQ_CORE_DIR}/SimDataEvents.cpp
    ${VQ_CORE_DIR}/SimDataEvents.h
    ${VQ_CORE_DIR}/Simulation.cpp
    ${VQ_CORE_DIR}/Simulation.h
    )

# Add executable called "vq" that is built from the source files
ADD_EXECUTABLE (vq main.cpp ${VQ_BINARY_DIR}/config.h ${VQ_MISC} ${VQ_SIMULATION} ${VQ_IO} ${VQ_CORE})
ADD_EXECUTABLE (mesher ${VQ_BINARY_DIR}/config.h mesher.cpp)

# Link the executables to the quakelib library.
TARGET_LINK_LIBRARIES (vq quakelib)
TARGET_LINK_LIBRARIES (mesher quakelib)

# Link the VQ simulator to the MPI library if found
IF(MPI_C_FOUND)
    TARGET_LINK_LIBRARIES (vq ${MPI_C_LIBRARIES})
    TARGET_LINK_LIBRARIES (vq ${MPI_CXX_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${MPI_C_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${MPI_CXX_LIBRARIES})
    SET_TARGET_PROPERTIES(vq PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(vq PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    SET_TARGET_PROPERTIES(mesher PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(mesher PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    message(STATUS "MPIEXEC=${MPIEXEC}")
ENDIF(MPI_C_FOUND)

IF(HDF5_FOUND)
    TARGET_LINK_LIBRARIES (vq ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
ENDIF(HDF5_FOUND)



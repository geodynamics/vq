# Make sure the compiler can find include files from the quakelib library. 
INCLUDE_DIRECTORIES(${VC_SOURCE_DIR}/quakelib/src)

# Set VC source directories
SET(VC_CORE_DIR ${VC_SOURCE_DIR}/src/core)
SET(VC_IO_DIR ${VC_SOURCE_DIR}/src/io)
SET(VC_MISC_DIR ${VC_SOURCE_DIR}/src/misc)
SET(VC_SIM_DIR ${VC_SOURCE_DIR}/src/simulation)

# and include files from VC itself
INCLUDE_DIRECTORIES (${VC_CORE_DIR})
INCLUDE_DIRECTORIES (${VC_IO_DIR})
INCLUDE_DIRECTORIES (${VC_MISC_DIR})
INCLUDE_DIRECTORIES (${VC_SIM_DIR})

# Make sure the linker can find the quakelib library once it is built. 
LINK_DIRECTORIES (${VC_BINARY_DIR}/quakelib)

# Define sets of the source files
SET(VC_MISC
    ${VC_MISC_DIR}/ConfigFile.cpp
    ${VC_MISC_DIR}/ConfigFile.h
    ${VC_MISC_DIR}/GreensFunctions.cpp
    ${VC_MISC_DIR}/GreensFunctions.h
    ${VC_MISC_DIR}/MPIDebugOutputStream.cpp
    ${VC_MISC_DIR}/MPIDebugOutputStream.h
    )

SET(VC_SIMULATION
    ${VC_SIM_DIR}/BASSAftershocks.cpp
    ${VC_SIM_DIR}/BASSAftershocks.h
    ${VC_SIM_DIR}/GracefulQuit.cpp
    ${VC_SIM_DIR}/GracefulQuit.h
    ${VC_SIM_DIR}/GreensInit.cpp
    ${VC_SIM_DIR}/GreensInit.h
    ${VC_SIM_DIR}/GreensKillInteraction.cpp
    ${VC_SIM_DIR}/GreensKillInteraction.h
    ${VC_SIM_DIR}/ProgressMonitor.cpp
    ${VC_SIM_DIR}/ProgressMonitor.h
    ${VC_SIM_DIR}/RunEvent.cpp
    ${VC_SIM_DIR}/RunEvent.h
    ${VC_SIM_DIR}/SanityCheck.cpp
    ${VC_SIM_DIR}/SanityCheck.h
    ${VC_SIM_DIR}/UpdateBlockStress.cpp
    ${VC_SIM_DIR}/UpdateBlockStress.h
    )

SET(VC_IO
    ${VC_IO_DIR}/CheckpointFileOutput.cpp
    ${VC_IO_DIR}/CheckpointFileOutput.h
    ${VC_IO_DIR}/CheckpointFileParse.cpp
    ${VC_IO_DIR}/CheckpointFileParse.h
    ${VC_IO_DIR}/EventOutput.cpp
    ${VC_IO_DIR}/EventOutput.h
    ${VC_IO_DIR}/GreensFileOutput.cpp
    ${VC_IO_DIR}/GreensFileOutput.h
    ${VC_IO_DIR}/HDF5Data.cpp
    ${VC_IO_DIR}/HDF5Data.h
    ${VC_IO_DIR}/ReadModelFile.cpp
    ${VC_IO_DIR}/ReadModelFile.h
    )

SET(VC_CORE
    ${VC_CORE_DIR}/SimError.h
    ${VC_CORE_DIR}/SimFramework.cpp
    ${VC_CORE_DIR}/SimFramework.h
    ${VC_CORE_DIR}/SimTimer.cpp
    ${VC_CORE_DIR}/SimTimer.h
    ${VC_CORE_DIR}/VCBlock.cpp
    ${VC_CORE_DIR}/VCBlock.h
    ${VC_CORE_DIR}/VCComm.cpp
    ${VC_CORE_DIR}/VCComm.h
    ${VC_CORE_DIR}/VCCommPartition.cpp
    ${VC_CORE_DIR}/VCCommPartition.h
    ${VC_CORE_DIR}/VCEvent.cpp
    ${VC_CORE_DIR}/VCEvent.h
    ${VC_CORE_DIR}/VCInitBlocks.cpp
    ${VC_CORE_DIR}/VCInitBlocks.h
    ${VC_CORE_DIR}/VCParams.cpp
    ${VC_CORE_DIR}/VCParams.h
    ${VC_CORE_DIR}/VCSimData.cpp
    ${VC_CORE_DIR}/VCSimData.h
    ${VC_CORE_DIR}/VCSimDataBlocks.cpp
    ${VC_CORE_DIR}/VCSimDataBlocks.h
    ${VC_CORE_DIR}/VCSimDataEvents.cpp
    ${VC_CORE_DIR}/VCSimDataEvents.h
    ${VC_CORE_DIR}/VCSimulation.cpp
    ${VC_CORE_DIR}/VCSimulation.h
    )

# Add executable called "vc" that is built from the source files
ADD_EXECUTABLE (vc main.cpp ${VC_MISC} ${VC_SIMULATION} ${VC_IO} ${VC_CORE})
ADD_EXECUTABLE (mesher mesher.cpp)

# Link the executables to the quakelib library.
TARGET_LINK_LIBRARIES (vc quakelib)
TARGET_LINK_LIBRARIES (mesher quakelib)

# Link the VC simulator to the MPI library if found
IF(MPI_C_FOUND)
    TARGET_LINK_LIBRARIES (vc ${MPI_C_LIBRARIES})
    TARGET_LINK_LIBRARIES (vc ${MPI_CXX_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${MPI_C_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${MPI_CXX_LIBRARIES})
    SET_TARGET_PROPERTIES(vc PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(vc PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    SET_TARGET_PROPERTIES(mesher PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(mesher PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
ENDIF(MPI_C_FOUND)

IF(HDF5_FOUND)
    TARGET_LINK_LIBRARIES (vc ${HDF5_LIBRARIES})
    TARGET_LINK_LIBRARIES (mesher ${HDF5_LIBRARIES})
ENDIF(HDF5_FOUND)



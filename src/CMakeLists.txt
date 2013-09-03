# Make sure the compiler can find include files from the quakelib library. 
INCLUDE_DIRECTORIES (${VC_SOURCE_DIR}/src/core)
INCLUDE_DIRECTORIES (${VC_SOURCE_DIR}/src/io)
INCLUDE_DIRECTORIES (${VC_SOURCE_DIR}/src/misc)
INCLUDE_DIRECTORIES (${VC_SOURCE_DIR}/src/simulation)

# Make sure the linker can find the Hello library once it is built. 
LINK_DIRECTORIES (${VC_BINARY_DIR}/quakelib)

# Define sets of the source files
SET(VC_MISC misc/ConfigFile.cpp misc/ConfigFile.h misc/GreensFunctions.cpp misc/GreensFunctions.h misc/MPIDebugOutputStream.cpp misc/MPIDebugOutputStream.h)
SET(VC_SIMULATION simulation/AddAsperities.cpp simulation/AddAsperities.h simulation/AddNoise.cpp simulation/AddNoise.h simulation/BadFaultKill.cpp simulation/BadFaultKill.h simulation/BASSAftershocks.cpp simulation/BASSAftershocks.h simulation/BGPoissonEvents.cpp simulation/BGPoissonEvents.h simulation/BlockValCompute.h simulation/BlockValCompute.cpp simulation/DepthDepVelocity.cpp simulation/DepthDepVelocity.h simulation/EventRecorder.cpp simulation/EventRecorder.h simulation/GracefulQuit.cpp simulation/GracefulQuit.h simulation/GreensInit.cpp simulation/GreensInit.h simulation/GreensKillInteraction.cpp simulation/GreensKillInteraction.h simulation/ProgressMonitor.cpp simulation/ProgressMonitor.h simulation/RunEvent.cpp simulation/RunEvent.h simulation/SanityCheck.cpp simulation/SanityCheck.h simulation/UpdateBlockStress.cpp simulation/UpdateBlockStress.h)
SET(VC_IO io/CheckpointFileOutput.cpp io/CheckpointFileOutput.h io/CheckpointFileParse.cpp io/CheckpointFileParse.h io/EqSimFileOutput.cpp io/EqSimFileOutput.h io/EqSimFileParse.cpp io/EqSimFileParse.h io/GreensFileOutput.cpp io/GreensFileOutput.h io/HDF5Data.cpp io/HDF5Data.h io/HDF5DataShare.cpp io/HDF5DataShare.h io/SystemFileOutput.cpp io/SystemFileOutput.h)
SET(VC_CORE core/SimError.h core/SimFramework.cpp core/SimFramework.h core/SimTimer.cpp core/SimTimer.h core/VCBlock.cpp core/VCBlock.h core/VCComm.cpp core/VCComm.h core/VCCommPartition.cpp core/VCCommPartition.h core/VCCommSpecExec.cpp core/VCCommSpecExec.h core/VCEvent.cpp core/VCEvent.h core/VCInitBlocks.cpp core/VCInitBlocks.h core/VCParams.cpp core/VCParams.h core/VCSimData.cpp core/VCSimData.h core/VCSimDataBlocks.cpp core/VCSimDataBlocks.h core/VCSimDataEvents.cpp core/VCSimDataEvents.h core/VCSimulation.cpp core/VCSimulation.h)

# Add executable called "vc" that is built from the source files
ADD_EXECUTABLE (vc main.cpp ${VC_MISC} ${VC_SIMULATION} ${VC_IO} ${VC_CORE})
#ADD_EXECUTABLE (mesher mesher.cpp)

# Link the executables to the quakelib library.
TARGET_LINK_LIBRARIES (vc quakelib)
#TARGET_LINK_LIBRARIES (mesher quakelib)

# Link the VC simulator to the MPI library if found
IF(DEFINED MPI_C_FOUND)
    TARGET_LINK_LIBRARIES (vc ${MPI_C_LIBRARIES})
    TARGET_LINK_LIBRARIES (vc ${MPI_CXX_LIBRARIES})
    SET_TARGET_PROPERTIES(vc PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(vc PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
ENDIF(DEFINED MPI_C_FOUND)

IF(DEFINED HDF5_FOUND)
    TARGET_LINK_LIBRARIES (vc ${HDF5_LIBRARIES})
ENDIF(DEFINED HDF5_FOUND)



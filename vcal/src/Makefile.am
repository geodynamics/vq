VC_MISC = \
	$(srcdir)/misc/ConfigFile.cpp \
	$(srcdir)/misc/ConfigFile.h \
	$(srcdir)/misc/GreensFunctions.cpp \
	$(srcdir)/misc/GreensFunctions.h \
	$(srcdir)/misc/MPIDebugOutputStream.cpp \
	$(srcdir)/misc/MPIDebugOutputStream.h

VC_SIMULATION = \
	$(srcdir)/simulation/AddAsperities.cpp \
	$(srcdir)/simulation/AddAsperities.h \
	$(srcdir)/simulation/AddNoise.cpp \
	$(srcdir)/simulation/AddNoise.h \
	$(srcdir)/simulation/BadFaultKill.cpp \
	$(srcdir)/simulation/BadFaultKill.h \
	$(srcdir)/simulation/BASSAftershocks.cpp \
	$(srcdir)/simulation/BASSAftershocks.h \
	$(srcdir)/simulation/BGPoissonEvents.cpp \
	$(srcdir)/simulation/BGPoissonEvents.h \
	$(srcdir)/simulation/BlockValCompute.cpp \
	$(srcdir)/simulation/BlockValCompute.h \
	$(srcdir)/simulation/DepthDepVelocity.cpp \
	$(srcdir)/simulation/DepthDepVelocity.h \
	$(srcdir)/simulation/EventRecorder.cpp \
	$(srcdir)/simulation/EventRecorder.h \
	$(srcdir)/simulation/GracefulQuit.cpp \
	$(srcdir)/simulation/GracefulQuit.h \
	$(srcdir)/simulation/GreensInit.cpp \
	$(srcdir)/simulation/GreensInit.h \
	$(srcdir)/simulation/GreensKillInteraction.cpp \
	$(srcdir)/simulation/GreensKillInteraction.h \
	$(srcdir)/simulation/ProgressMonitor.cpp \
	$(srcdir)/simulation/ProgressMonitor.h \
	$(srcdir)/simulation/RunEvent.cpp \
	$(srcdir)/simulation/RunEvent.h \
	$(srcdir)/simulation/SanityCheck.cpp \
	$(srcdir)/simulation/SanityCheck.h \
	$(srcdir)/simulation/UpdateBlockStress.cpp \
	$(srcdir)/simulation/UpdateBlockStress.h

VC_IO = \
	$(srcdir)/io/CheckpointFileOutput.cpp \
	$(srcdir)/io/CheckpointFileOutput.h \
	$(srcdir)/io/CheckpointFileParse.cpp \
	$(srcdir)/io/CheckpointFileParse.h \
	$(srcdir)/io/EqSimFileOutput.cpp \
	$(srcdir)/io/EqSimFileOutput.h \
	$(srcdir)/io/EqSimFileParse.cpp \
	$(srcdir)/io/EqSimFileParse.h \
	$(srcdir)/io/GreensFileOutput.cpp \
	$(srcdir)/io/GreensFileOutput.h \
	$(srcdir)/io/HDF5Data.cpp \
	$(srcdir)/io/HDF5Data.h \
	$(srcdir)/io/HDF5DataShare.cpp \
	$(srcdir)/io/HDF5DataShare.h \
	$(srcdir)/io/SystemFileOutput.cpp \
	$(srcdir)/io/SystemFileOutput.h

VC_CORE = \
	$(srcdir)/core/SimError.h \
	$(srcdir)/core/SimFramework.cpp \
	$(srcdir)/core/SimFramework.h \
	$(srcdir)/core/SimTimer.cpp \
	$(srcdir)/core/SimTimer.h \
	$(srcdir)/core/VCBlock.cpp \
	$(srcdir)/core/VCBlock.h \
	$(srcdir)/core/VCComm.cpp \
	$(srcdir)/core/VCComm.h \
	$(srcdir)/core/VCCommPartition.cpp \
	$(srcdir)/core/VCCommPartition.h \
	$(srcdir)/core/VCCommSpecExec.cpp \
	$(srcdir)/core/VCCommSpecExec.h \
	$(srcdir)/core/VCEvent.cpp \
	$(srcdir)/core/VCEvent.h \
	$(srcdir)/core/VCInitBlocks.cpp \
	$(srcdir)/core/VCInitBlocks.h \
	$(srcdir)/core/VCParams.cpp \
	$(srcdir)/core/VCParams.h \
	$(srcdir)/core/VCSimData.cpp \
	$(srcdir)/core/VCSimData.h \
	$(srcdir)/core/VCSimDataBlocks.cpp \
	$(srcdir)/core/VCSimDataBlocks.h \
	$(srcdir)/core/VCSimDataEvents.cpp \
	$(srcdir)/core/VCSimDataEvents.h \
	$(srcdir)/core/VCSimulation.cpp \
	$(srcdir)/core/VCSimulation.h

bin_PROGRAMS = vc
vc_CXXFLAGS = $(CXXFLAGS) $(OPENMP_CFLAGS) $(HDF5_CPPFLAGS) -I$(srcdir)/../quakelib/src -I$(srcdir)/core -I$(srcdir)/io -I$(srcdir)/simulation -I$(srcdir)/graphics -I$(srcdir)/misc
vc_LDFLAGS = $(HDF5_LDFLAGS) -L$(srcdir)/../quakelib/.libs -lquake-1.0
vc_SOURCES = main.cpp $(VC_MISC) $(VC_SIMULATION) $(VC_IO) $(VC_CORE)


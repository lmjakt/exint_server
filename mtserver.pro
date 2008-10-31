TEMPLATE = app
CONFIG = qt release pgsql thread
HEADERS = 	version.h \
		server/server.h \
		server/connectionobject.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		stat/stat.h \
		netArray/netArray.h \
		server/processor.h \
		server/anovaProcessor.h \
		server/euclidSortProcessor.h \
		server/kClusterProcess.h \
		raw/dataStructs.h \
		server/experimentCompareProcess.h \
		server/flatExptCompare.h \
		protocol/protocol.h \
		experiment/experiment.h \
		server/blastClient.h \
		util/dataExtractor.h \
		util/normaliser.h \
		util/ProbeStats.h \
		util/pathTracer.h \
		util/experimentTracer.h \
		util/sorted_floats.h 
SOURCES = 	main.cpp \
		server/server.cpp \
		server/connectionobject.cpp \
		raw/probe_set.cpp \
		raw/probeSetSet2.cpp \
		stat/stat.cpp \
		netArray/netArray.cpp \
		server/processor.cpp \
		server/anovaProcessor.cpp \
		server/euclidSortProcessor.cpp \
		server/kClusterProcess.cpp \
		raw/dataStructs.cpp \
		server/experimentCompareProcess.cpp \
		server/flatExptCompare.cpp \
		protocol/protocol.cpp \
		experiment/experiment.cpp \
		server/blastClient.cpp \
		util/dataExtractor.cpp \
		util/normaliser.cpp \
		util/ProbeStats.cpp \
		util/pathTracer.cpp \
		util/experimentTracer.cpp
TARGET = mtserver_2.12


#############################################################################
# Makefile for building: mtserver_2.12
# Generated by qmake (1.07a) (Qt 3.3.3) on: Fri Apr  1 15:01:21 2005
# Project:  mtserver.pro
# Template: app
# Command: $(QMAKE) -o Makefile mtserver.pro
#############################################################################

####### Compiler, tools and options

CC       = gcc
CXX      = g++
LEX      = flex
YACC     = yacc
CFLAGS   = -pipe -O2 -g -pipe -m64 -D_REENTRANT  -DQT_NO_DEBUG -DQT_THREAD_SUPPORT
CXXFLAGS = -pipe -O2 -g -pipe -m64 -D_REENTRANT  -DQT_NO_DEBUG -DQT_THREAD_SUPPORT -DHAVE_NAMESPACE_STD -DHAVE_CXX_STRING_HEADER -DDLLIMPORT=""
LEXFLAGS = 
YACCFLAGS= -d
INCPATH  = -I/usr/lib64/qt-3.3/mkspecs/default -I. -I$(QTDIR)/include -I/usr/local/pgsql/include
LINK     = g++
LFLAGS   = 
LIBS     = $(SUBLIBS) -L$(QTDIR)/lib -L/usr/X11R6/lib64  -L/usr/local/pgsql/lib -lqt-mt -lXext -lX11 -lm -lpthread -lpq++
AR       = ar cqs
RANLIB   = 
MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic
QMAKE    = qmake
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p

####### Output directory

OBJECTS_DIR = ./

####### Files

HEADERS = server/server.h \
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
SOURCES = main.cpp \
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
OBJECTS = main.o \
		server.o \
		connectionobject.o \
		probe_set.o \
		probeSetSet2.o \
		stat.o \
		netArray.o \
		processor.o \
		anovaProcessor.o \
		euclidSortProcessor.o \
		kClusterProcess.o \
		dataStructs.o \
		experimentCompareProcess.o \
		flatExptCompare.o \
		protocol.o \
		experiment.o \
		blastClient.o \
		dataExtractor.o \
		normaliser.o \
		ProbeStats.o \
		pathTracer.o \
		experimentTracer.o
FORMS = 
UICDECLS = 
UICIMPLS = 
SRCMOC   = server/moc_server.cpp
OBJMOC = moc_server.o
DIST	   = mtserver.pro
QMAKE_TARGET = mtserver_2.12
DESTDIR  = 
TARGET   = mtserver_2.12

first: all
####### Implicit rules

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(UICDECLS) $(OBJECTS) $(OBJMOC)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJMOC) $(OBJCOMP) $(LIBS)

mocables: $(SRCMOC)
uicables: $(UICDECLS) $(UICIMPLS)

$(MOC): 
	( cd $(QTDIR)/src/moc && $(MAKE) )

Makefile: mtserver.pro  /usr/lib64/qt-3.3/mkspecs/default/qmake.conf 
	$(QMAKE) -o Makefile mtserver.pro
qmake: 
	@$(QMAKE) -o Makefile mtserver.pro

dist: 
	@mkdir -p .tmp/mtserver_2.12 && $(COPY_FILE) --parents $(SOURCES) $(HEADERS) $(FORMS) $(DIST) .tmp/mtserver_2.12/ && ( cd `dirname .tmp/mtserver_2.12` && $(TAR) mtserver_2.12.tar mtserver_2.12 && $(GZIP) mtserver_2.12.tar ) && $(MOVE) `dirname .tmp/mtserver_2.12`/mtserver_2.12.tar.gz . && $(DEL_FILE) -r .tmp/mtserver_2.12

mocclean:
	-$(DEL_FILE) $(OBJMOC)
	-$(DEL_FILE) $(SRCMOC)

uiclean:

yaccclean:
lexclean:
clean: mocclean
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) $(TARGET)


FORCE:

####### Compile

main.o: main.cpp server/server.h \
		server/connectionobject.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h \
		raw/dataStructs.h

server.o: server/server.cpp server/server.h \
		server/connectionobject.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o server.o server/server.cpp

connectionobject.o: server/connectionobject.cpp server/connectionobject.h \
		server/anovaProcessor.h \
		server/euclidSortProcessor.h \
		server/kClusterProcess.h \
		server/experimentCompareProcess.h \
		server/flatExptCompare.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		raw/dataStructs.h \
		netArray/netArray.h \
		protocol/protocol.h \
		experiment/experiment.h \
		util/dataExtractor.h \
		util/normaliser.h \
		util/ProbeStats.h \
		util/experimentTracer.h \
		server/blastClient.h \
		util/pathTracer.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o connectionobject.o server/connectionobject.cpp

probe_set.o: raw/probe_set.cpp raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o probe_set.o raw/probe_set.cpp

probeSetSet2.o: raw/probeSetSet2.cpp raw/probeSetSet2.h \
		server/server.h \
		raw/probe_set.h \
		raw/dataStructs.h \
		server/connectionobject.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o probeSetSet2.o raw/probeSetSet2.cpp

stat.o: stat/stat.cpp stat/stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o stat.o stat/stat.cpp

netArray.o: netArray/netArray.cpp netArray/netArray.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o netArray.o netArray/netArray.cpp

processor.o: server/processor.cpp server/processor.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o processor.o server/processor.cpp

anovaProcessor.o: server/anovaProcessor.cpp server/anovaProcessor.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		server/connectionobject.h \
		raw/dataStructs.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o anovaProcessor.o server/anovaProcessor.cpp

euclidSortProcessor.o: server/euclidSortProcessor.cpp server/euclidSortProcessor.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		server/connectionobject.h \
		raw/dataStructs.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o euclidSortProcessor.o server/euclidSortProcessor.cpp

kClusterProcess.o: server/kClusterProcess.cpp server/kClusterProcess.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		server/connectionobject.h \
		raw/dataStructs.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/blastClient.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o kClusterProcess.o server/kClusterProcess.cpp

dataStructs.o: raw/dataStructs.cpp raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dataStructs.o raw/dataStructs.cpp

experimentCompareProcess.o: server/experimentCompareProcess.cpp server/experimentCompareProcess.h \
		raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o experimentCompareProcess.o server/experimentCompareProcess.cpp

flatExptCompare.o: server/flatExptCompare.cpp server/flatExptCompare.h \
		raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o flatExptCompare.o server/flatExptCompare.cpp

protocol.o: protocol/protocol.cpp protocol/protocol.h \
		netArray/netArray.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o protocol.o protocol/protocol.cpp

experiment.o: experiment/experiment.cpp experiment/experiment.h \
		netArray/netArray.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o experiment.o experiment/experiment.cpp

blastClient.o: server/blastClient.cpp server/blastClient.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o blastClient.o server/blastClient.cpp

dataExtractor.o: util/dataExtractor.cpp util/dataExtractor.h \
		util/normaliser.h \
		raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dataExtractor.o util/dataExtractor.cpp

normaliser.o: util/normaliser.cpp util/normaliser.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o normaliser.o util/normaliser.cpp

ProbeStats.o: util/ProbeStats.cpp util/ProbeStats.h \
		util/dataExtractor.h \
		util/normaliser.h \
		raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ProbeStats.o util/ProbeStats.cpp

pathTracer.o: util/pathTracer.cpp util/pathTracer.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pathTracer.o util/pathTracer.cpp

experimentTracer.o: util/experimentTracer.cpp util/experimentTracer.h \
		util/pathTracer.h \
		util/dataExtractor.h \
		raw/probe_set.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o experimentTracer.o util/experimentTracer.cpp

moc_server.o: server/moc_server.cpp  server/server.h server/connectionobject.h \
		raw/probe_set.h \
		raw/probeSetSet2.h \
		netArray/netArray.h \
		util/pathTracer.h \
		server/kClusterProcess.h \
		server/blastClient.h \
		raw/dataStructs.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_server.o server/moc_server.cpp

server/moc_server.cpp: $(MOC) server/server.h
	$(MOC) server/server.h -o server/moc_server.cpp

####### Install

install:  

uninstall:  

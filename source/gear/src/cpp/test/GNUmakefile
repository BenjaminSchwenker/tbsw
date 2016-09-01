######################################################
#
# GNUmakefile to build the GEAR library examples
#
# @author F. Gaede, DESY
# @version $Id: GNUmakefile,v 1.7 2007-05-23 10:36:21 gaede Exp $
######################################################


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar


INSTALL = $(GEARBASE)/lib
BINDIR = $(GEARBASE)/bin
INCLUDE =  $(GEARINCLUDE) 


CPPFLAGS = $(GEARCPPFLAGS) $(INCLUDE) 

LIBS = -L$(INSTALL) -L$(SIOINSTALL)/lib  

ifdef GEAR_USE_CGA
include $(GEARBASE)/src/cpp/gearcga/src/GNUmakefile.cga
LIBS += -lgearcga $(CGALIBS)
CPPFLAGS += -DCGA
endif

LIBS += -lgearxml -lgear

#----- testong tpcproto ---------
ifdef GEAR_USE_AIDA
  CPPFLAGS += -D GEAR_USE_AIDA `aida-config --include`
  LIBS += `aida-config --lib`
  LIBS += `root-config --glibs`
endif
#--------------------------------


BINS := $(patsubst %.cc,$(BINDIR)/%,$(wildcard *.cc))


bin: $(BINS)

clean:
	rm -f $(BINS)



$(BINDIR)/testgear: testgear.cc $(INSTALL)/libgear.a
	$(CXX) -o $(BINDIR)/testgear $(CPPFLAGS) testgear.cc  $(LIBS) 

$(BINDIR)/mergeXML: mergeXML.cc $(INSTALL)/libgear.a
	$(CXX) -o $(BINDIR)/mergeXML $(CPPFLAGS) mergeXML.cc  $(LIBS)




$(BINDIR)/testtpcproto: testtpcproto.cc $(INSTALL)/libgear.a
	$(CXX) -o $(BINDIR)/testtpcproto $(CPPFLAGS) testtpcproto.cc  $(LIBS) 

$(BINDIR)/testVTXgear: testVTXgear.cc $(INSTALL)/libgear.a
	$(CXX) -o $(BINDIR)/testVTXgear $(CPPFLAGS) testVTXgear.cc  $(LIBS) 

$(BINDIR)/testSiPlanesgear: testSiPlanesgear.cc $(INSTALL)/libgear.a
	$(CXX) -o $(BINDIR)/testSiPlanesgear $(CPPFLAGS) testSiPlanesgear.cc  $(LIBS) 

######################################################
#
# GNUmakefile to build the LCCD library test examples
#
# @author F. Gaede, DESY
# @version $Id: GNUmakefile,v 1.6 2006-11-14 17:34:27 gaede Exp $
######################################################
SIOINSTALL = $(LCCDBASE)/sio

# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar


INSTALL = $(LCCDBASE)/lib
BINDIR = $(LCCDBASE)/test

INCLUDE =  $(LCCDINCLUDE)
CPPFLAGS = $(LCCDCPPFLAGS) $(LCCDINCLUDE) 
CPPFLAGS += -I $(LCIO)/src/cpp/include

LIBS = -L $(INSTALL) -llccd -L$(LCIO)/lib -llcio -L$(LCIO)/sio/lib -lsio -lz


BINSRC = readcalfile.cc createcalfile.cc testdbfile.cc readdatafile.cc multiblefolders.cc

ifdef CondDBMySQL
BINSRC += calfiletodb.cc calfilefromdb.cc tagdbfolder.cc createdbfile.cc
LIBS += -L $(CondDBMySQL)/lib -lconddb
endif

BINS := $(patsubst %.cc,$(BINDIR)/%,$(BINSRC))


objects := $(patsubst %.cc,%.o,$(filter-out $(BINSRC), $(wildcard *.cc))  )


# --- check dependencies (from S.Aplin) 
dependence = $(patsubst %.cc,%.d,$(wildcard *.cc))

%.d: %.cc
	set -e; $(CC) -MM $(CPPFLAGS) $< | sed -e 's/\($*\)\.o[ :]*/\1.o $*.d : /g' > $@; \
	[ -s $@ ] || rm -f $@

include $(dependence)


bin: $(BINS)

clean:
	rm -f $(BINS) *.o *.d


$(BINDIR)/readcalfile: readcalfile.cc $(INSTALL)/liblccd.a $(objects)
	mkdir -p $(BINDIR)
	$(CXX) -o $(BINDIR)/readcalfile $(CPPFLAGS) readcalfile.cc  $(objects) $(LIBS) 

$(BINDIR)/createcalfile: createcalfile.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/createcalfile $(CPPFLAGS) createcalfile.cc   $(objects)  $(LIBS) 

$(BINDIR)/testdbfile: testdbfile.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/testdbfile $(CPPFLAGS) testdbfile.cc   $(objects)  $(LIBS) 

$(BINDIR)/readdatafile: readdatafile.cc $(INSTALL)/liblccd.a $(objects)
	mkdir -p $(BINDIR)
	$(CXX) -o $(BINDIR)/readdatafile $(CPPFLAGS) readdatafile.cc  $(objects) $(LIBS) 



ifdef CondDBMySQL
$(BINDIR)/calfiletodb: calfiletodb.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/calfiletodb $(CPPFLAGS) calfiletodb.cc   $(objects)  $(LIBS) 

$(BINDIR)/multiblefolders: multiblefolders.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/multiblefolders $(CPPFLAGS) multiblefolders.cc   $(objects)  $(LIBS) 

$(BINDIR)/calfilefromdb: calfilefromdb.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/calfilefromdb $(CPPFLAGS) calfilefromdb.cc   $(objects)  $(LIBS) 

$(BINDIR)/tagdbfolder: tagdbfolder.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/tagdbfolder $(CPPFLAGS) tagdbfolder.cc   $(objects)  $(LIBS) 

$(BINDIR)/createdbfile: createdbfile.cc $(INSTALL)/liblccd.a $(objects)
	$(CXX) -o $(BINDIR)/createdbfile $(CPPFLAGS) createdbfile.cc   $(objects)  $(LIBS) 
endif
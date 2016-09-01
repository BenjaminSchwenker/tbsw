######################################################
#
# GNUmakefile to build the GEAR library
#
# @author Frank Gaede
# @version $Id: GNUmakefile,v 1.1.1.1 2005-09-13 14:41:03 gaede Exp $
######################################################


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar

INSTALL  = $(GEARBASE)/lib
INCLUDE  = $(GEARINCLUDE)
CPPFLAGS = $(GEARCPPFLAGS) $(INCLUDE)


objects := $(patsubst %.cc,%.o,$(wildcard *.cc))


$(INSTALL)/libgearxml.a: $(objects)
	$(AR) cr $(INSTALL)/libgearxml.a  $?


lib: $(INSTALL)/libgearxml.a


clean:
	rm -f *.o $(INSTALL)/libgearxml.a





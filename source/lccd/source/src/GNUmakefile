######################################################
#
# GNUmakefile to build the LCCD library
#
# @author Frank Gaede, DESY
# @version $Id: GNUmakefile,v 1.2 2005-02-11 11:05:27 gaede Exp $
######################################################


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar


INSTALL =  $(LCCDBASE)/lib
INCLUDE =  $(LCCDINCLUDE)
CPPFLAGS = $(LCCDCPPFLAGS) $(INCLUDE)


# --- check dependencies (from S.Aplin) 
dependence = $(patsubst %.cc,%.d,$(wildcard *.cc))

%.d: %.cc
	set -e; $(CC) -MM $(CPPFLAGS) $< | sed -e 's/\($*\)\.o[ :]*/\1.o $*.d : /g' > $@; \
	[ -s $@ ] || rm -f $@

include $(dependence)



objects := $(patsubst %.cc,%.o,$(wildcard *.cc))

$(INSTALL)/liblccd.a: $(objects)
	mkdir -p $(INSTALL)
	$(AR) cr $(INSTALL)/liblccd.a  $?

#lib: $(objects)
lib: $(INSTALL)/liblccd.a

clean:
	rm -f *.o *.d





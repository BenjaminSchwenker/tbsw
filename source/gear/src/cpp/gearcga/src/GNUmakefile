##
## $Id: GNUmakefile,v 1.1 2006-05-29 13:03:18 musat Exp $
## $Name: not supported by cvs2svn $
##
## History
## first implementation of Gear Distance and Point Properties with the 
## Mokka Common Geometry Access (CGA) by 
## Gabriel Musat (musat@poly.in2p3.fr), May 2006
##
## see CGA documentation at 
## http://polype.in2p3.fr/geant4/tesla/www/mokka/
##        software/doc/CGADoc/CGAIndex.html
##-------------------------------------------------------

# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif
                                                                                
AR = ar

CPPFLAGS += \
	-I $(MOKKA)/source/Geometry/CGA/include \
	$(GEARCPPFLAGS) $(GEARINCLUDE)

INSTALL  = $(GEARBASE)/lib

objects := $(patsubst %.cc,%.o,$(wildcard *.cc))
                                                                                
$(INSTALL)/libgearcga.a: $(objects)
	$(AR) cr $(INSTALL)/libgearcga.a  $? 

lib: $(INSTALL)/libgearcga.a
                                                                                
clean:
	rm -f *.o $(INSTALL)/libgearcga.a


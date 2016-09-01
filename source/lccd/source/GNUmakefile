######################################################
#
# GNUmakefile to build the LCCD library
#
# @author F. Gaede, DESY
# @version $Id: GNUmakefile,v 1.4 2005-03-11 13:54:18 gaede Exp $
######################################################

LCCDBASE := $(PWD)/..
export LCCDBASE 


#CLHEP_INCLUDE := $(CLHEP)/include


LCCDINCLUDE := -I $(LCCDBASE)/source/include -I $(LCIO)/src/cpp/include
#-I $(CLHEP_INCLUDE)
export LCCDINCLUDE

ifdef LCCDDEBUG
 LCCDCPPFLAGS = -g -Wall -ansi
else
 LCCDCPPFLAGS = -O3 -Wall -ansi
endif


ifdef CondDBMySQL
 LCCDCPPFLAGS += -D USE_CONDDB -I $(CondDBMySQL)/include
endif

# to set the same default as in ConditionsDBMySQL, do
# strings lib/libconddb.so | grep -E ".+:.+:.+:.+"
ifdef LCCD_DB_INIT
  LCCDCPPFLAGS += -D LCCD_DB_INIT_DEFAULT=\"$(LCCD_DB_INIT)\"
else
  LCCDCPPFLAGS += -D LCCD_DB_INIT_DEFAULT=\"localhost:lccd_test:calvin:hobbes\"
endif

export LCCDCPPFLAGS


.PHONY: lib clean CORE doc

all: lib bin


lib: CORE


bin: lib
	$(MAKE) -C ./test bin

CORE:
	 echo $(LCCDBASE) ; $(MAKE) -C ./src lib


clean:
	$(MAKE) -C ./src clean; $(MAKE) -C ./test clean; \
	rm -rf $(LCCDBASE)/lib/liblccd.a

doc:
	cd $(LCCDBASE)/doc ; doxygen





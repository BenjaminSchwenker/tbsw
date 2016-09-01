###################################################################
#
# GNUmakefile to build the gear library
#
# @author F. Gaede, DESY
# @version $Id: GNUmakefile,v 1.5 2006-05-29 13:03:18 musat Exp $
###################################################################


GEARBASE := $(GEAR)
export GEARBASE 

GEARINCLUDE := -I $(GEARBASE)/src/cpp/include
export GEARINCLUDE

ifdef GEARDEBUG
 GEARCPPFLAGS = -g -pg -Wall -pedantic -Wno-long-long -ansi 
else
 GEARCPPFLAGS = -O3 -Wall -pedantic -Wno-long-long -ansi 
endif


export GEARCPPFLAGS

packages := gearimpl gearxml

ifdef GEAR_USE_CGA
packages += gearcga
endif

subdirs := $(patsubst %, $(GEARBASE)/src/cpp/%, $(packages) )

includedirs := $(patsubst %, $(GEARBASE)/src/cpp/include/%, $(packages) )


.PHONY: lib clean bin doc


all: lib bin


$(includedirs):
	@echo $(includedirs) ; \
	cd $(GEARBASE)/src/cpp/include ; \
	ln -fs $(patsubst $(GEARBASE)/src/cpp/include/%, ../%/include, $@ ) $(patsubst $(GEARBASE)/src/cpp/include/%, %, $@ )

#	ln -s $(patsubst $(GEARBASE)/src/cpp/include/%, $(GEARBASE)/src/cpp/%/include, $@ ) $@




lib: $(includedirs)
	@for i in $(packages); do \
	if [ -f "$(GEARBASE)/src/cpp/$$i/src/GNUmakefile" ] ; then  \
	echo "Building library for $$i..."; \
	(cd $(GEARBASE)/src/cpp/$$i/src; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) lib); fi ; \
	done

#	echo "ln -s $(GEARBASE)/src/cpp/$$i/include  $(GEARBASE)/src/cpp/include/$$i" ; \
#	ln -s $(GEARBASE)/src/cpp/$$i/include  $(GEARBASE)/src/cpp/include/$$i ; \


bin: lib $(GEARBASE)/lib/libgear.a  $(GEARBASE)/lib/libgearxml.a
	mkdir -p $(GEARBASE)/bin ; \
	$(MAKE) -C ./test bin


clean:
	@for i in $(subdirs); do \
	if [ -f "$$i/src/GNUmakefile" ] ; then \
	echo "cleaning up in $$i..."; \
	(cd $$i/src; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) clean); fi ; done ; \
	rm $(GEARBASE)/bin/* ; \
	rm -rf $(GEARBASE)/src/cpp/include/*.h $(includedirs)


doc:
	cd $(GEARBASE)/doc ; doxygen







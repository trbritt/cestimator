SUBDIRS := src
BINDIR := bin 

all: $(SUBDIRS) 
$(SUBDIRS):
	mkdir -p $(BINDIR)
	$(MAKE) -C $@
	ln -sf $@/cestimator ./bin/cestimator

.PHONY: all $(SUBDIRS) clean

clean: 
	$(MAKE) -C $(SUBDIRS) clean
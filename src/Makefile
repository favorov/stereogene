makefiles=StereoGene.mak \
	Projector.mak \
	Confounder.mak \
	Binning.mak \
	Smoother.mak

objdir=./obj

.PHONY:	all clean

all: $(objdir) ccvars $(makefiles)
	$(MAKE) -f StereoGene.mak
	$(MAKE) -f Projector.mak 
	$(MAKE) -f Confounder.mak
	$(MAKE) -f Binning.mak 
	$(MAKE) -f Smoother.mak
	
$(objdir):
	mkdir -p $(objdir)

clean:
	$(MAKE) -f StereoGene.mak clean
	$(MAKE) -f Projector.mak clean
	$(MAKE) -f Confounder.mak clean
	$(MAKE) -f Binning.mak clean 
	$(MAKE) -f Smoother.mak clean


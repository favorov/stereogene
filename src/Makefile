#****************************************************************************#
#SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001
#$Id: Makefile 1900 2013-07-09 07:36:01Z favorov $
#****************************************************************************#
#this file can be started only from parent directory

makefiles=\
	StereoGene.mak \
	Smoother.mak \
	Projector.mak \
	Binner.mak \
	Confounder.mak

objdir=./obj

.PHONY:	all clean

all: $(objdir) ccvars $(makefiles)
	$(MAKE) -f StereoGene.mak
	$(MAKE) -f Smoother.mak
	$(MAKE) -f Projector.mak 
	$(MAKE) -f Binner.mak 
	$(MAKE) -f Confounder.mak
	$(MAKE) -f ParseGenes.mak
	
$(objdir):
	mkdir -p $(objdir)

clean:
	$(MAKE) -f StereoGene.mak clean
	$(MAKE) -f Smoother.mak clean
	$(MAKE) -f Projector.mak clean
	$(MAKE) -f Binner.mak clean
	$(MAKE) -f Confounder.mak clean
	$(MAKE) -f ParseGenes.mak clean


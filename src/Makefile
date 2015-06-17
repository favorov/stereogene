#****************************************************************************#
#SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001
#$Id: Makefile 1900 2013-07-09 07:36:01Z favorov $
#****************************************************************************#
#this file can be started only from parent directory

makefile=StereoGene.mak
objdir=./obj

.PHONY:	all clean zip cleanzip


all: $(objdir) ccvars $(makefile) 
	$(MAKE) -f $(makefile)

$(objdir):
	mkdir -p $(objdir)

clean:
	$(MAKE) -f $(makefile) clean

zip:
	$(MAKE) -f $(makefile) zip

cleanzip:
	$(MAKE) -f $(makefile) cleanzip


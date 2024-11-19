#use gmake!
name=Stereogene
#bibname=gene-best-friends
#bibname=$(name)

tex="xelatex"

.PHONY: all clean

all: $(name).pdf

$(name).pdf: $(name).tex 
	#$(name).bbl	
	$(tex) $(name).tex
	$(tex) $(name).tex

#$(name).bbl: $(bibname).bib
#	pdflatex $(name).tex
#	bibtex $(name)	

clean: 
	rm -f $(name).bbl $(name).blg $(name).log $(name).aux $(name).dvi *.gz *~



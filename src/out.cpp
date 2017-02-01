/*
 * out.cpp
 *
 *  Created on: 12.12.2013
 *      Author: 1
 */
#include "track_util.h"
#include <sys/file.h>

statTest *MannW;
int XYCorrScale=100;
//================================================================= write results in ENCODE Broad Peak format
void printCorrelations(){
	verb("\nWrite correlations...\n");
	char b[1024];

	errStatus="printCorrelations";
	//============================================================= Write the foreground distribution
	if(writeDistr) printFgDistr();
	if(writeBPeak) printBroadPeak();
	if(writeDistCorr)
		printChrDistances(strcat(strcpy(b,outFile),".dist"));
	if(outSpectr) XYfgCorrelation.printSpect(strcat(strcpy(b,outFile),".spect"));
	if(outChrom ) printChomosomes(strcat(strcpy(b,outFile),".chrom"));
	errStatus=0;
}

void printChrDistances(char *fname){
	FILE *f=xopen(fname,"wt");
	XYbgcorrelation.normilize();
	XYfgCorrelation.normilize();
	normChromDist();
	fprintf(f,"# %s vs %s\n",bTrack1.name, bTrack2.name);

	//========================================
	double *autoCorrX=bTrack1.autoCorr, *autoCorrY=bTrack2.autoCorr;
	fprintf(f,"x");
	if(doAutoCorr) fprintf(f,"\tAutoCorr1\tAutoCorr2");
	fprintf(f,"\tBkg\tFg\tFgPlus\tFgMinus");

	//========================================
	if(outChrom){
		for(int i=0; i<n_chrom; i++){
			Chromosome *chr=chrom_list+i;
			if(chr->densCount)	fprintf(f,"\t%s",chr->chrom);
		}
	}
	fprintf(f,"\n");
	for(int j=0; j<profWithFlanksLength; j++){
		int k=(j+profWithFlanksLength/2)%profWithFlanksLength;
		fprintf(f,"%i",(j-profWithFlanksLength/2)*binSize);
		if(doAutoCorr) fprintf(f,"\t%9.5f\t%9.5f",autoCorrX[k],autoCorrY[k]);
		fprintf(f,"\t%9.5f\t%9.5f\t%9.5f\t%9.5f", XYbgcorrelation.correlation[k]*XYCorrScale,
				XYfgCorrelation.correlation[k]*XYCorrScale, XYfgCorrelation.corrPlus[k]*XYCorrScale, XYfgCorrelation.corrMinus[k]*XYCorrScale);
		if(outChrom){
			for(int i=0; i<n_chrom; i++){
				Chromosome *chr=chrom_list+i;
				if(chr->densCount){
					if(chr->count > 1){
						double x=chr->distDens[k];
						fprintf(f,"\t%7.3f",x*XYCorrScale);
					}
					else{
						fprintf(f, "\tNA");
					}
				}
			}
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

void printChomosomes(char *fname){
	FILE *f=xopen(fname,"wt");
	fprintf(f,"%s\n%s\n",bTrack1.name, bTrack2.name);
	fprintf(f,"%6s\t%5s \t%5s \t%5s \t%5s\n",
			"chrom","av1", "av2", "cc",  "count");
	for(int i=0; i<n_chrom; i++){
		Chromosome *chr=chrom_list+i;
		if(chr->count > 1)
			fprintf(f,"%6s\t%6.2f\t%6.2f\t%6.2f\t%6.0f\n",
				chr->chrom, chr->av1/chr->count, chr->av2/chr->count,
				chr->corr/chr->count, chr->count);
	}

	fclose(f);
}
//======================================================== write report to cummulative report file
void getStat(double *set, int n, double &av, double &sd){
	av=sd=0;
	for(int i=0; i<n; i++){
		av+=set[i]; sd+=set[i]*set[i];
	}
	av/=n; sd=sqrt((sd-av*av*n)/(n-1));
}

double avBg=0,sdBg=0,avFg=0,sdFg=0;

void printStat(){
	verb("Write statistics\n");
	writeLog("Write statistics\n");
	char b[2048];
	FILE *f=0;

	MannW=MannWhitney(FgSet, nFg, BkgSet, nBkg);	// do Mann-Whitney test
	xverb("p-val=%e\nnWindows=%i\n=================================\n",
		MannW->pVal, nFg);

	getStat(FgSet, nFg, avFg,sdFg);
	getStat(BkgSet,nBkg,avBg,sdBg);
	const char *pc= pcorProfile==0 ? "NONE" : pcorProfile;

	bool fg=fileExists(statFileName);
	if((outRes & TAB)!=0) {
		f=gopen(statFileName,"a+t");
		if(f!=0) flockFile(f);
		else {
			fprintf(stderr,"Can not open file %s Error code=%i\n",statFileName, errno);
			writeLog("Can not open file %s\n",statFileName);
		}
	}
	if(!fg && f){								//================ write the header
		fprintf(f,"%-6s\t%-20s%-20s\t%-20s","id","version","name1","name2");
		fprintf(f,"\t%-6s\t%-6s\t%-6s\t%-6s","window","kern","nFgr","nBkg");
		fprintf(f,"\t%-8s\t%-8s\t\%-8s","Fg_Corr","Fg_av_Corr","FgCorr_sd");
		fprintf(f,"\t%-8s\t%-8s\t\%-8s","Bg_Corr","Bg_av_Corr","BgCorr_sd");
		fprintf(f,"\t%-8s\t%-7s", "Mann-Z","p-value");
		fprintf(f,"\t%-6s\n", "P-corr");
	}
	//==================================================== write the statistics
	if(f){
		char nm1[1024], nm2[1024];
		alTable.convert(bTrack1.name, nm1);
		alTable.convert(bTrack2.name, nm2);
		fprintf(f,"%08lx\t%-10s\t%-10s\t%-10s",id,version, nm1, nm2);
		fprintf(f,"\t%-6i\t%-6s\t%6i\t%6i",wSize,getKernelType(),nFg, nBkg);
		fprintf(f,"\t%8.4f\t%8.4f\t%8.4f",totCorr, avFg,sdFg);
		fprintf(f,"\t%8.4f\t%8.4f\t%8.4f",BgTotal, avBg, sdBg);
		fprintf(f,"\t%8.4f\t%-7.2e\t%s\n",MannW->z,MannW->pVal, pc);
		funlockFile(f);
		fclose(f);
	}

	//================================================== write parameters
	fg=fileExists(paramsFileName);
	if((outRes&TAB)!=0){
		f=gopen(paramsFileName,"a+t");
		if(f) flockFile(f);
		else {
			fprintf(stderr,"Can not open file %s  Error code=%i\n",paramsFileName,errno);
			writeLog("Can not open file %s\n",paramsFileName);
		}
	}
	if(!fg && f){								//================ write the header
		printParamNames(f);
		fprintf(f,"\n");
	}
	//==================================================== write the parameters
	if(f){
		printParams(f);
		fprintf(f,"\n");
		fclose(f);
	}
	if((outRes & XML)!=0) {
		sprintf(b,"%s.xml",statFileName);
		fg=fileExists(b);
		FILE *xml=0;
		if(!fg) {xml=xopen(b,"w"); fprintf(xml,"<xml>\n");}
		else{
			xml=xopen(b,"r+");
			fseek(xml,-7,SEEK_END);
		}
		flockFile(xml);
		fprintf(xml,"<run id=\"%08lx\" ver=\"%s\">\n", id, version);
		fprintf(xml,"\t<input track1=\"%s\" track2=\"%s\"/>\n",bTrack1.name,bTrack2.name);
		fprintf(xml,"\t<output out=\"%s\"/>\n",outFile);
		fprintf(xml,"\t<prm ");
		printXMLparams(xml);
		fprintf(xml,"/>\n");

		fprintf(xml,"\t<res ");
		fprintf(xml,"nFg=\"%i\" ",nFg);
		fprintf(xml,"nBkg=\"%i\" ", nBkg);
		fprintf(xml,"totCorr=\"%.4f\" ", totCorr);
		fprintf(xml,"FgAvCorr=\"%.4f\" ", FgAvCorr);
		fprintf(xml,"BgAvCorr=\"%.4f\" ", BgAvCorr);
		fprintf(xml,"MannZ=\"%.4f\" ",MannW->z);
		fprintf(xml,"pVal=\"%.2e\" ",MannW->pVal);

		fprintf(xml,"sdBg=\"%.4f\" ",sdBg);
		fprintf(xml,"sdFg=\"%.4f\" ",sdFg);
		fprintf(xml,"/>\n");
		fprintf(xml,"</run>\n");
		fprintf(xml,"</xml>\n");
		funlockFile(xml);
		fclose(xml);
	}
	writeLog("Write Statistics -> Done\n");
}

void XYCorrelation::printSpect(char *fname){
	FILE *f=xopen(fname,"wt");
	fprintf(f,"#\t%s\t%s\n",bTrack1.name, bTrack2.name);
	fprintf(f,"Wave_Length\tSpectrum1\tSpectrum2\n");
	double dx=0,dy=0;
	for(int i=0; i<profWithFlanksLength; i++){
		dx+=spectrumX[i]; dy+=spectrumY[i];
	}
	for(int i=0; i<profWithFlanksLength; i++){
		spectrumX[i]=sqrt(spectrumX[i]/dx*profWithFlanksLength);
		spectrumY[i]=sqrt(spectrumY[i]/dy*profWithFlanksLength);
	}

	for(int i=0; i<profWithFlanksLength/2; i++){
		double l=(double)(profWithFlanksLength)*binSize/(i+1);
		fprintf(f,"%9.2f\t%9.4f\t%9.4f\n",l,spectrumX[i],spectrumY[i]);
	}

	fclose(f);
}

//============================================== write the foreground distribution of the correlations
void printFgDistr(){
	char b[1024];
	FILE *fFgDistr=0;
	strcat(strcpy(b,outFile),".fg");
	fFgDistr=xopen(b,"wt");
	ScoredRange gp;

	for(int i=0; i<nPairs; i++){
		filePos2Pos(pairs[i].profPos,&gp,wSize);
		fprintf(fFgDistr,"%s\t%ld\t%ld\t%f\n",gp.chrom, gp.beg,gp.end, pairs[i].d);			// write the distribution: correlation, p-value, q-value
	}
	fclose(fFgDistr);
}

//============================================== write the BroadPeak file
void printBroadPeak(){
	char b[1024];
	ScoredRange *pPair;
	getMem(pPair,(nPairs+10), "printBroadPeak");
	int nPpint=0;
	//============================================================= Select pairs to print
	for(int i=0; i<nPairs; i++){
		PairEntry *pe=pairs+i;
		double pv= bgHist.pValm(pe->d);						// get p--value
		double qv=pv/(fgHist.pValp(pairs[i].d));			//q-value=nPairs*pv/(nPairs-pe->rank);
		if(wSize > wStep) qv=qv*wSize/wStep;				// correct q-value if step < window
		if(qv > 0.99999) qv=0.99999;
		pv=-log10(pv); qv=-log10(qv);						// log p-value & q-value
		if(pv > pVal && qv > qVal){							// p-value and q-value are good enough
			filePos2Pos(pe->profPos,&pPair[nPpint],wSize);		// get chromosome position
			pPair[nPpint].score=pe->d;						// define information for write
			nPpint++;
		}
	}


	//============================================================== Clear overlaps
	for(int i=1; i<nPpint; i++){
		if(strcmp(pPair[i-1].chrom,pPair[i].chrom)==0){
			if(pPair[i-1].end > pPair[i].beg){			//===== pairs overpep
				if(pPair[i-1].beg > pPair[i].beg)
					pPair[i].beg=pPair[i-1].end;
				else if(pPair[i-1].score < pPair[i].score)	//===== select highest score
					pPair[i-1].end=pPair[i].beg;
				else
					pPair[i].beg=pPair[i-1].end;
			}
		}
	}

	//=============================================================== Print broad peak file
	strcat(strcpy(b,outFile),".bpeak");
	FILE *fBpeak=xopen(b,"wt");
	//============================================================= write header
	fprintf(fBpeak,"track type=broadPeak name=%s__%s\n",bTrack1.name,bTrack2.name);
	for(int i=0; i<nPpint; i++){
		ScoredRange *pe=pPair+i;
		double pv= bgHist.pValm(pe->score);			// get p-value
		double qv=pv/(fgHist.pValp(pairs[i].d));  	//qv=nPairs*pv/(nPairs-pe->rank);// get q-value
		if(wSize > wStep) qv=qv*wSize/wStep;
		if(qv > 0.99999) qv=0.99999;
		if(pv==0) pv=1.e-256;
		if(qv==0) qv=1.e-256;
		pv=-log10(pv); qv=-log10(qv);							// log p-value & q-value
		int score=int((pe->score+1)*500);
		if(pe->end - pe->beg > 2)								// not-empty interval
		fprintf(fBpeak,"%s\t%li\t%li\t.\t%i\t.\t.\t%5.1f\t%5.1f\n",
				pe->chrom, pe->beg, pe->end, score,pv,qv);
	}
	xfree(pPair,"pPair");
	fclose(fBpeak);
}
void printRmd(){
	char b[2048];
	sprintf(b,"%sreport_r_template.Rmd",resPath);
	
	if(!fileExists(b))
	{
		FILE *f=xopen(b,"wt");
		fprintf(f, "---	\n");
		fprintf(f, "title: \"Report\"	\n");
		fprintf(f, "output: html_document	\n");
		fprintf(f, "params:	\n");
		fprintf(f, " track1: !r as.character(\"\")	\n");
		fprintf(f, " track2: !r as.character(\"\")	\n");
		fprintf(f, " pc: !r as.character(\"\")	\n");
		fprintf(f, " name: !r as.character(\"\")	\n");
		fprintf(f, " window: !r NA	\n");
		fprintf(f, " kernel: !r NA	\n");
		fprintf(f, " nFgr: !r NA	\n");
		fprintf(f, " nBkg: !r NA	\n");
		fprintf(f, " Bkg_av: !r NA	\n");
		fprintf(f, " Fg_av: !r NA	\n");
		fprintf(f, " Bkg_sd: !r NA	\n");
		fprintf(f, " Fg_sd: !r NA	\n");
		fprintf(f, " tot_cor: !r NA	\n");
		fprintf(f, " avCorr: !r NA	\n");
		fprintf(f, " Mann_Z: !r NA	\n");
		fprintf(f, " p_value: !r NA	\n");
		fprintf(f, "\n");
		fprintf(f, "---	\n");
		fprintf(f, "```{r echo=FALSE}	\n");
		fprintf(f, "window <- params$window	\n");
		fprintf(f, "kernel <- params$kernel	\n");
		fprintf(f, "nFgr <- params$nFgr	\n");
		fprintf(f, "nBkg <- params$nBkg	\n");
		fprintf(f, "Bkg_av <- params$Bkg_av	\n");
		fprintf(f, "Fg_av <- params$Fg_av	\n");
		fprintf(f, "Bkg_sd <- params$Bkg_sd   	\n");
		fprintf(f, "Fg_sd <- params$Fg_sd	\n");
		fprintf(f, "tot_cor <- params$tot_cor	\n");
		fprintf(f, "avCorr <- params$avCorr	\n");
		fprintf(f, "Mann_Z <- params$Mann_Z  	\n");
		fprintf(f, "p_value <- params$p_value	\n");
		fprintf(f, "```	\n");
		fprintf(f, "\n");
		fprintf(f, "```{r eval=(params$track1!=\"\"), echo=FALSE, comment=\"\", results=\'asis\'}	\n");
		fprintf(f, "	cat(paste(\"<p word-break: break-all>track1: \", params$track1, \"</p>\",\"<p>track2: \", params$track2, \"</p>\", sep=\"\"))	\n");
		fprintf(f, "```	\n");
		fprintf(f, "```{r eval=(params$pc!=\"\"), echo=FALSE, comment=\"\", results=\'asis\'}	\n");
		fprintf(f, "	cat(paste(\"<p word-break: break-all>partial correlation track: \", params$pc, \"</p>\", sep=\"\"))	\n");
		fprintf(f, "```	\n");
		fprintf(f, "\n");
		fprintf(f, "\n");
		fprintf(f, "Parameter | Value  	\n");
		fprintf(f, "------------- | -------------  	\n");
		fprintf(f, "window | `r window`  	\n");
		fprintf(f, "kernel | `r kernel` 	\n");
		fprintf(f, "nFgr | `r nFgr`	\n");
		fprintf(f, "Fg_av | `r Fg_av`	\n");
		fprintf(f, "Fg_sd | `r Fg_sd`	\n");
		fprintf(f, "nBkg | `r nBkg`	\n");
		fprintf(f, "Bkg_av | `r Bkg_av`	\n");
		fprintf(f, "Bkg_sd | `r Bkg_sd`	\n");
		fprintf(f, "tot_cor | `r tot_cor`	\n");
		fprintf(f, "avCorr | `r avCorr`	\n");
		fprintf(f, "Mann_Z | `r Mann_Z`	\n");
		fprintf(f, "p_value | `r p_value`	\n");
		fprintf(f, "\n");
		fprintf(f, "```{r eval=(params$name!=\"\"), echo=FALSE}	\n");
		fprintf(f, "name <- params$name	\n");
		fprintf(f, "\n");
		fprintf(f, "fg <- read.table(paste(name, \'.fg\', sep = \'\'))	\n");
		fprintf(f, "bkg<- read.table(paste(name, \'.bkg\', sep = \'\'))	\n");
		fprintf(f, "dist <- read.table(paste(name, \'.dist\', sep = \'\'), header=TRUE)	\n");
		fprintf(f, "#  Define plot limits	\n");
		fprintf(f, "\n");
		fprintf(f, "y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y))	\n");
		fprintf(f, "y_lim2 <- c(min(min(dist$Fg),min(dist$Fg)),max(max(dist$Fg),max(dist$Fg)))	\n");
		fprintf(f, "\n");
		fprintf(f, "x_lim2 <- c(-10000,10000)	\n");
		fprintf(f, "\n");
		fprintf(f, "# set x scale to kilobases	\n");
		fprintf(f, "x_lim2 <- x_lim2/1000	\n");
		fprintf(f, "\n");
		fprintf(f, "#get chromosome data for plots, example for chr1. 	\n");
		fprintf(f, "#Some times you should also reset y_lim for plots	\n");
		fprintf(f, "#fg_chrom <- fg[fg[,1]==\"chr1\",]	\n");
		fprintf(f, "#dist_chrom <- dist$chr1	\n");
		fprintf(f, "\n");
		fprintf(f, "\n");
		fprintf(f, "# save plot to pdf	\n");
		fprintf(f, "#  create the plot	\n");
		fprintf(f, "old.par <- par( no.readonly = TRUE )	\n");
		fprintf(f, "par( mfrow = c( 2, 1 ), oma = c( 0, 0, 0, 0 ),mar=c(3,3,2,1),mgp=c(1.6,0.45,0))	\n");
		fprintf(f, "\n");
		fprintf(f, "\n");
		const char*dens="density";
		if(XYCorrScale!=1) sprintf(b,"%s*%i",dens,XYCorrScale); else strcpy(b,dens);
		fprintf(f, "plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), xlab=\'correlation coefficient\',ylab=\'%s\',	\n",b);
		fprintf(f, "col=\'red\', main=\'Distribution of correlations\',	\n");
		fprintf(f, "cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2)	\n");
		fprintf(f, "lines(density(fg[,4]), col=\'blue\', lwd=2)	\n");
		fprintf(f, "\n");
		fprintf(f, "#plot line for chomosome	\n");
		fprintf(f, "#lines(density(fg[,4]), col=\'green\', lwd=2)	\n");
		fprintf(f, "\n");
		fprintf(f, "\n");
		fprintf(f, "plot(dist$x/1000, dist$Fg, type=\'l\',col=\'blue\', ylim=y_lim2, xlim=x_lim2,	\n");
		fprintf(f, "main=\'Cross-correlation function\',xlab=\'Distance (kb)\',ylab=\'density*100\',cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2)	\n");
		fprintf(f, "lines(dist$x/1000,dist$Bkg , col=\'red\',lwd=2)	\n");
		fprintf(f, "#plot line for chomosome	\n");
		fprintf(f, "#lines(dist$x/1000, dist_chrom , col=\'green\',lwd=2)	\n");
		fprintf(f, "\n");
		fprintf(f, "par( old.par )	\n");
		fprintf(f, "```	\n");

		fclose(f);		
	}
	
}
void printRreport(){
	char *s,b[2048], fname[1024], fn[1024];
	alTable.convert(outFile,fn);

	if(sdFg==0) getStat(FgSet,nFg,avFg,sdFg);
	if(sdBg==0) getStat(BkgSet,nBkg,avBg,sdBg);

	strcat(strcpy(b,fn),"_report.r");
	FILE *f=xopen(b,"wt");

	strcpy(b,outFile);
	s=strrchr(b,'/'); if(s==0) s=outFile; else s++; strcpy(fname,s);

	
	fprintf(f, "library(\"markdown\")\n");
	
	fprintf(f, "args = commandArgs(TRUE)\n");
	fprintf(f, "fname1<-\"%s\"\n", bTrack1.name);
	fprintf(f, "fname2<-\"%s\"\n", bTrack2.name);
	if (pcorProfile!=0) {
		fprintf(f,"pc_fname <- \"%s\"\n",pcorProfile);
	}
	else {
		fprintf(f,"pc_fname <- \"\"\n");
	}
	fprintf(f, "\nif (length(args)>=2) {\n");
  	fprintf(f, "#track names\n");
  	fprintf(f, "	fname1 <- args[1]\n"); 
  	fprintf(f, "	fname2 <- args[2]\n");
	fprintf(f, "} \n\n");
	fprintf(f, "if (length(args)==3){\n");
  	fprintf(f, "#partial correlation track\n");
  	fprintf(f, "	pc_fname <- args[3]\n");
	fprintf(f, "} \n\n");

	fprintf(f, "rmarkdown::render(\"report_r_template.Rmd\", \"html_document\", \n");
    fprintf(f, "              params=list(\n");
  	fprintf(f, "track1=fname1, \n");
  	fprintf(f, "track2=fname2, \n");
  	fprintf(f, "pc=pc_fname, \n");
  	fprintf(f, "name=\"%s\", \n", fname);
  	fprintf(f, "window=\"%i\", \n", wSize);
  	fprintf(f, "kernel=\"%s\",\n", getKernelType());
  	fprintf(f, "nFgr=\"%i\",\n", nFg);
  	fprintf(f, "nBkg=\"%i\",\n", nBkg);
  	fprintf(f, "Bkg_av=\"%.4f\",\n", BgAvCorr);
  	fprintf(f, "Fg_av=\"%.4f\",\n", FgAvCorr);
  	fprintf(f, "Bkg_sd=\"%.4f\", \n", sdBg);
  	fprintf(f, "Fg_sd=\"%.4f\",\n", sdFg);
  	fprintf(f, "tot_cor=\"%.4f\",\n", totCorr);
//  	fprintf(f, "avCorr=\"%.4f\",\n", FgAvCorr);
  	fprintf(f, "Mann_Z=\"%.4f\",  \n", MannW->z);
  	fprintf(f, "p_value=\"%.2e\" \n", MannW->pVal);	
  	fprintf(f, "), output_file = file.path(getwd(), \"%s.html\"))\n", fname);
	
	fclose(f);
}

void printR(){
	char fn[1024], *s, b[2048], fname[1024];
	alTable.convert(outFile, fn);
//	const char *cex="cex.axis = 0.8,  cex.lab = 0.8,  cex.main = 0.8", *lwd="lwd=2",
//			*lab="xlab=\'correlation coefficient\',ylab=\'density\'";
//
//	int    x0,x1;
//	double y0,y1;
//	correlation.getLimits(x0,x1, y0, y1);
//
	strcat(strcpy(b,fn),".r");
	FILE *f=xopen(b,"wt");

	strcpy(b,outFile);
	s=strrchr(b,'/'); if(s==0) s=outFile; else s++; strcpy(fname,s);
//==============================================================================
	fprintf(f," #  Read the data  \n");
	fprintf(f," name <-  \'%s\'  \n\n",fname);
	fprintf(f," fg <- read.table(paste(name, '.fg', sep = '')) \n");
	fprintf(f," bkg<- read.table(paste(name, '.bkg', sep = '')) \n");
	fprintf(f," dist <- read.table(paste(name, '.dist', sep = ''), header=TRUE) \n\n");
	fprintf(f," #  Define plot limits \n\n");
	fprintf(f," y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y)) \n");
	fprintf(f," x_lim2 <- c(-%i,%i) \n\n",crossWidth,crossWidth);
	fprintf(f," # set x scale to kilobases \n");
	fprintf(f," x_lim2 <- x_lim2/1000 \n\n");
	fprintf(f," #get chromosome data for plots, example for chr1.  \n");
	fprintf(f," #Some times you should also reset y_lim for plots \n");
	fprintf(f," #fg_chrom <- fg[fg[,1]==\"chr1\",] \n");
	fprintf(f," #dist_chrom <- dist$chr1 \n\n\n");
	fprintf(f," # save plot to pdf \n");
	if(writePDF) fprintf(f," pdf(paste(name,'.pdf', sep=''), height = 9, width = 5) \n\n");
	fprintf(f," #  create the plot \n");
	fprintf(f," old.par <- par( no.readonly = TRUE ) \n");
	fprintf(f," par( mfrow = c( %i, 1 ), oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0)) \n\n",
			outWIG?3:2);
	char sub[1024]; if(writePDF) sub[0]=0; else sprintf(sub,"\\n%s",fname);
	fprintf(f," plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), xlab='correlation coefficient',ylab='density', \n");
	fprintf(f," col='red', main='Distribution of correlations%s', \n",sub);
	fprintf(f," cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2) \n");
	fprintf(f," lines(density(fg[,4]), col='blue', lwd=2) \n");
	fprintf(f," #plot line for chomosome \n");
	fprintf(f," #lines(density(fg[,4]), col='green', lwd=2) \n\n\n");
	fprintf(f," plot(dist$x/1000, dist$Fg, type='l',col='blue', xlim=x_lim2, \n");
	const char*dens="density";
	if(XYCorrScale!=1) sprintf(b,"%s*%i",dens,XYCorrScale); else strcpy(b,dens);
	fprintf(f," main='Cross-correlation function',xlab='Distance (kb)',ylab='%s',cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2) \n",b);
	fprintf(f," lines(dist$x/1000,dist$Bkg , col='red',lwd=2) \n");
	fprintf(f," #plot line for chomosome \n");
	fprintf(f," #lines(dist$x/1000, dist_chrom , col='green',lwd=2) \n\n");
	if(outWIG){
		const char *ss="";
		if(LCScale==LOG_SCALE) ss="(log)";
		if(LCScale==LOG_LOG_SCALE) ss="(log log)";
		double nn=fiveFDR*2; if(nn > 1000) nn=1000;
		fprintf(f,"lc=read.table(paste(name, \'.LChist\', sep = \'\'),, header=TRUE)\n");
		fprintf(f,"xmax=%.0f; x0=xmax-150; x1=x0+20; y0=90; dy=7;\n",nn);
		fprintf(f,"plot(lc$val, lc$cdf_obs*100000, col='blue', lwd=2, \n");
		fprintf(f,"type=\'l\',ylim=c(0,100), xlim=c(0,xmax), ylab=\'FDR\',");
		fprintf(f,"xlab=\'normalized LC-values %s\',\n",ss);
		fprintf(f,"main=\'Local correlation distributions\');\n");
		fprintf(f,"lines(lc$val, lc$cdf_exp*100000, col=\'red\', lwd=2);\n");
		fprintf(f,"lines(lc$val, lc$FDR, col='black', lwd=3);\n");
		fprintf(f," text(x0, y0+7,name,pos=4); y0=y0-5\n");
		fprintf(f,"segments(x0,y0,x1,y0,col='red',lwd=2); text(x1+5,y0,'expected (1 - cdf) * 100000',pos=4)\n");
		fprintf(f,"y0=y0-dy;\n");
		fprintf(f,"segments(x0,y0,x1,y0,col='blue',lwd=2); text(x1+5,y0,'observed (1 - cdf) * 100000',pos=4)\n");
		fprintf(f,"y0=y0-dy;\n");
		fprintf(f,"segments(x0,y0,x1,y0,col='black',lwd=3); text(x1+5,y0,'FDR',pos=4)\n");
		fprintf(f,"y0=y0-dy;\n");
		fprintf(f,"segments(x0,y0,x1,y0,col='gray',lwd=3); text(x1+5,y0,'FDR=5%%',pos=4)\n");

		fprintf(f,"segments(0,5,xmax,5,col='gray',lwd=1); text(x1+5,y0,'FDR',pos=4)\n");
	}
	fprintf(f," par( old.par ) \n\n");
	if(writePDF) fprintf(f," dev.off() \n");

//===========================   OLD Version ==============================
//	fprintf(f,"#  Read the data \n\n");
//
//	fprintf(f,"fg <- read.table(\'%s.fg\')\n",fname);
//	fprintf(f,"bkg<- read.table(\'%s.bkg\')\n",fname);
//
//	fprintf(f,"dist <- read.table(\'%s.dist\', header=TRUE)\n",fname);
//	fprintf(f,"\n#  Define plot limits\n\n");
//
//	fprintf(f,"y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y))\n");
//	fprintf(f,"y_lim2 <- c(%.2f,%.2f)\n",y0,y1);
//
//	fprintf(f,"x_lim2 <- c(%i,%i)\n\n",x0,x1);
//
//	fprintf(f,"# set x scale to kilobases\n");
//	fprintf(f,"x_lim2 <- x_lim2/1000\n");
//
//	fprintf(f,"\n#  create the plot\n\n");
//
//	fprintf(f,"old.par <- par( no.readonly = TRUE )\n");
//	fprintf(f,"par( mfrow = c( 2, 1 ), oma = c( 0.5, 0, 2, 0 ),mar=c(2.5,3,1.5,1),mgp=c(1,0.3,0))\n\n");
//
//	fprintf(f,"plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), %s,\n",lab);
//	fprintf(f,"col=\'red\', main=\'Distribution of correlations, p-value=%.1e\',\n",MannW->pVal );
//	fprintf(f,"%s,%s)\n\n",cex,lwd);
//	fprintf(f,"lines(density(fg[,4]), col=\'blue\', %s)\n\n",lwd);
//
//	fprintf(f,"plot(dist$x/1000, dist$Fg, type=\'l\',col=\'blue\', ylim=y_lim2, xlim=x_lim2,\n");
//	fprintf(f,"main=\'Cross-correlation function\',xlab=\'Distance (kb)\',ylab=\'density,%%\',%s,%s)\n",cex,lwd);
//	fprintf(f,"#lines(dist$x/1000,dist$FgPlus, col=\'cyan\',%s)\n",lwd);
//	fprintf(f,"#lines(dist$x/1000,dist$FgMinus, col=\'brown\',%s)\n",lwd);
//	fprintf(f,"lines(dist$x/1000,dist$Bkg , col=\'red\',%s)\n",lwd);
//
//	char fn1[1024], fn2[1024];
//	strcpy(fn1,fname);
//	s=strchr(fn1,'~');
//	if(s){*s=0; s++; strcpy(fn2,s);}
//	else {*fn2=0;}
//
//	fprintf(f,"title(\'%s\\n%s\',cex.main = 0.8, outer = TRUE)\n\n", fn1,fn2);
//	fprintf(f,"par( old.par )\n");

	fclose(f);
}


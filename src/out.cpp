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
	if(writeDistr) {
		printFgDistr();
		printBgDistr();
	}
	if(writeDistCorr)
		printChrDistances(strcat(strcpy(b,outFile),".dist"));
	if(outSpectr) XYfgCorrelation.printSpect(strcat(strcpy(b,outFile),".spect"));
	if(outChrom ) printChomosomes(strcat(strcpy(b,outFile),".chrom"));
	if(doAutoCorr) printAuto();
	errStatus=0;
}

void printAuto(){
	track1->writeAuto();
	track2->writeAuto();
}

void printChrDistances(char *fname){
	FILE *f=xopen(fname,"wt");
	XYbgcorrelation.normilize();
	XYfgCorrelation.normilize();
	normChromDist();
	fprintf(f,"# %s vs %s\n",track1->name, track2->name);

	//========================================
	fprintf(f,"x");
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
		int k=(j + profWithFlanksLength / 2 ) % profWithFlanksLength;
		fprintf(f,"%i",(j - profWithFlanksLength / 2) * binSize);
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
	fprintf(f,"%s\n%s\n",track1->name, track2->name);
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
	if(n>1) {
		sd=sqrt((sd-av*av/n)/(n-1));
	}
	else sd=FNA;
	if(n > 0) av/=n;
	else av=FNA;

}

double avBg=FNA,sdBg=FNA,avFg=FNA,sdFg=FNA;
void printStat(){
	verb("Write statistics\n");
	writeLog("Write statistics\n");
	char b[2048];
	MannW=MannWhitney(FgSet, nFg, BkgSet, nBkg);	// do Mann-Whitney test
	if(MannW ==0) {
		mannW_Z=FNA;
		mannW_p=FNA;
		xverb("p-val=NA\nnWindows=%i\n=================================\n",
				mannW_p, nFg);
	}
	else{
		mannW_Z=MannW->z;
		mannW_p=MannW->pVal;
		xverb("p-val=%e\nnWindows=%i\n=================================\n",
				mannW_p, nFg);
	}
	getStat(FgSet, nFg, avFg,sdFg);
	getStat(BkgSet,nBkg,avBg,sdBg);
	FILE *f=0;
	bool fg=true;

	if(statFileName){
		fg=fileExists(statFileName);
		if((outRes & TAB)!=0) {
			f=gopen(statFileName,"a+t");
			if(f!=0) flockFile(f);
			else {
				writeLogErr("Can not open file %s\n",statFileName);
			}
		}
		if(!fg && f){	//================ write the header
			printStatHeader(f);
		}
		//==================================================== write the statistics
		if(f){
			printStat(f);
			fclose(f);
		}
	}

	//================================================== write parameters
	if(paramsFileName){
		fg=fileExists(paramsFileName);
		if((outRes&TAB)!=0){
			f=gopen(paramsFileName,"a+t");
			if(f) flockFile(f);
			else {
				writeLogErr("Can not open file %s\n",paramsFileName);
			}
		}
		if(!fg && f){								//================ write the header
			printParamNames(f);
		}
		//==================================================== write the parameters
		if(f){
			printParams(f);
			funlockFile(f);
			fclose(f);
		}
	}
	if(statFileName){
		if((outRes & XML)!=0) {
			sprintf(b,"%s.xml",statFileName);
			fg=fileExists(b);
			FILE *xml=0;
			if(!fg) {xml=xopen(b,"wb"); fprintf(xml,"<xml>\n");}
			else{
				xml=xopen(b,"r+b");
				fseek(xml,-7,SEEK_END);
			}
			flockFile(xml);
			printXML(xml);
			fprintf(xml,"</xml>\n");
			funlockFile(xml);
			fclose(xml);
		}
	}
	if(customKern){
		FILE*cust=gopen("kernels","a+");
		fprintf(cust,"%s\t\"%s\"\n", printId(), customKern);
		fclose(cust);
	}
	writeLog("Write Statistics -> Done\n");
}

void XYCorrelation::printSpect(char *fname){
	FILE *f=xopen(fname,"wt");
	fprintf(f,"#\t%s\t%s\n",track1->name, track2->name);
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

//============================================== write the foreground and background distributions of the correlations
void printBgDistr(){
	char b[1024];
	strcat(strcpy(b,outFile),".bkg");					// open file for background observations
	FILE* fbkg=xopen(b,"wt");
	for(int i=0; i<nBkg; i++) fprintf(fbkg,"%f\n",BkgSet[i]);
	fclose(fbkg);
}

void printFgDistr(){
	char b[1024];
	FILE *fFgDistr=0;
	strcat(strcpy(b,outFile),".fg");
	fFgDistr=xopen(b,"w");
	ScoredRange gp;
	if(writeDistr==DISTR_DETAIL){
		for(int i=0; i<nFgPos; i++){
			filePos2Pos(FgCorr[i].profPos,&gp,wSize);
			fprintf(fFgDistr,"%s\t%ld\t%ld\t%f\n",gp.chrom, gp.beg,gp.end, FgCorr[i].d);			// write the distribution: correlation, p-value, q-value
		}
	}
	else{
		for(int i=0; i<nFg; i++) fprintf(fFgDistr,"%f\n",FgSet[i]);
	}
	fclose(fFgDistr);
}

//=====================================================
//=====================================================
//=====================================================
//=====================================================
//=====================================================
//=====================================================
const char *template1="report_r_template1.Rmd";
const char *template2="report_r_template2.Rmd";
const char *template3="report_r_template3.Rmd";

void printRmd(){

	char b[2048];
	int nPlot=1;
	
	if(writeDistr==DISTR_SHORT)  sprintf(b,"%s%s",resPath,template1);
	else if(!doAutoCorr)		{sprintf(b,"%s%s",resPath,template2); nPlot=2;}
	else						{sprintf(b,"%s%s",resPath,template3); nPlot=3;}

//	if(!fileExists(b))
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
		if(doAutoCorr){
			fprintf(f, "auto1 <- read.table(paste(fname1, \'.auto\', sep = \'\'), header=FALSE)	\n");
			fprintf(f, "auto2 <- read.table(paste(fname2, \'.auto\', sep = \'\'), header=FALSE)	\n");
			fprintf(f, "y_lim3 <- c(min(min(auto1$V2),min(auto2$V2)),max(max(auto1$V2),max(auto2$V2)))	\n");
		}
		fprintf(f, "#  Define plot limits	\n");
		fprintf(f, "\n");
		if(writeDistr==DISTR_SHORT){
			fprintf(f, "y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,1])$y))	\n");
		}
		else
			fprintf(f, "y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y))	\n");
		fprintf(f, "y_lim2 <- c(min(min(dist$Fg),min(dist$Fg)),max(max(dist$Fg),max(dist$Fg)))	\n");
		fprintf(f, "\n");
		fprintf(f, "x_lim2 <- c(-10000,10000)	\n");
		fprintf(f, "\n");
		fprintf(f, "# set x scale to kilobases	\n");
		fprintf(f, "x_lim2 <- x_lim2/1000	\n");
		fprintf(f, "\n");
		if(writeDistr==DISTR_DETAIL){
			fprintf(f, "#get chromosome data for plots, example for chr1. 	\n");
			fprintf(f, "#Some times you should also reset y_lim for plots	\n");
			fprintf(f, "#fg_chrom <- fg[fg[,1]==\"chr1\",]	\n");
			fprintf(f, "#dist_chrom <- dist$chr1	\n");
			fprintf(f, "\n");
			fprintf(f, "\n");
		}
		fprintf(f, "# save plot to pdf	\n");
		fprintf(f, "#  create the plot	\n");
		fprintf(f, "old.par <- par( no.readonly = TRUE )	\n");
		fprintf(f, "par( mfrow = c( %i, 1 ), oma = c( 0, 0, 0, 0 ),mar=c(3,3,2,1),mgp=c(1.6,0.45,0))\n",nPlot);
		fprintf(f, "\n");
		fprintf(f, "\n");
		const char*dens="density";
		if(XYCorrScale!=1) sprintf(b,"%s*%i",dens,XYCorrScale); else strcpy(b,dens);
		fprintf(f, "plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), xlab=\'correlation coefficient\',ylab=\'%s\',	\n",b);
		fprintf(f, "col=\'red\', main=\'Distribution of correlations\',	\n");
		fprintf(f, "cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2)	\n");
		if(writeDistr==DISTR_SHORT)
			fprintf(f, "lines(density(fg[,1]), col=\'blue\', lwd=2)	\n");
		else
			fprintf(f, "lines(density(fg[,4]), col=\'blue\', lwd=2)	\n");
		fprintf(f, "\n");
		fprintf(f, "#plot line for chomosome	\n");
		if(writeDistr==DISTR_SHORT)
			fprintf(f, "#lines(density(fg[,1]), col=\'green\', lwd=2)	\n");
		else
			fprintf(f, "#lines(density(fg[,4]), col=\'green\', lwd=2)	\n");
		fprintf(f, "\n\n");
		fprintf(f, "plot(dist$x/1000, dist$Fg, type=\'l\',col=\'blue\', ylim=y_lim2, xlim=x_lim2,	\n");
		fprintf(f, "main=\'Cross-correlation function\',xlab=\'Distance (kb)\',ylab=\'density*100\',cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2)	\n");
		fprintf(f, "lines(dist$x/1000,dist$Bkg , col=\'red\',lwd=2)	\n");
		if(doAutoCorr){
			fprintf(f, "x_lim3=c(max(min(auto1$V1/1000),-100), min(max(auto1$V1/1000),100))\n");
			fprintf(f, "plot( auto1$V1/1000, auto1$V2, type=\'l\',col=\'blue\',ylim=y_lim3, xlim=x_lim3, main=\'Autocorrelations\'");
					fprintf(f, ",xlab=\'Distance (kb)\',ylab=\'autocorr\',lwd=2)\n");
			fprintf(f, "lines(auto2$V1/1000, auto2$V2, type=\'l\',col=\'red\',lwd=2)\n");
			fprintf(f, "dy=y_lim3[2]-y_lim3[1];  y1=y_lim3[1]+dy*0.5;  y2=y_lim3[1]+dy*0.75\n");

			fprintf(f, "text(0,y1,fname1,col='blue',pos = 4)\n");
			fprintf(f, "text(0,y2,fname2,col='red',pos = 4)\n");
		}

		fprintf(f, "\n");
		fprintf(f, "par( old.par )	\n");
		fprintf(f, "```	\n");

		fclose(f);		
	}
	
}
void printRreport(){
	char *s,b[2048], fname[1024];

	if(sdFg==0) getStat(FgSet,nFg,avFg,sdFg);
	if(sdBg==0) getStat(BkgSet,nBkg,avBg,sdBg);

	strcat(strcpy(b,outFile),"_report.r");
	FILE *f=xopen(b,"wt");

	strcpy(b,outFile);
	s=strrchr(b,'/'); if(s==0) s=outFile; else s++; strcpy(fname,s);

	
	fprintf(f, "library(\"markdown\")\n");
	
	fprintf(f, "args = commandArgs(TRUE)\n");
	fprintf(f, "fname1<-\"%s\"\n", track1->name);
	fprintf(f, "fname2<-\"%s\"\n", track2->name);
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

	const char *xtemplate;
	if(doAutoCorr) 		 				xtemplate=template3;
	else if(writeDistr == DISTR_SHORT) 	xtemplate=template1;
	else				 				xtemplate=template2;

	fprintf(f, "rmarkdown::render(\"%s\", \"html_document\", \n",xtemplate);
    fprintf(f, "              params=list(\n");
  	fprintf(f, "track1=fname1, \n");
  	fprintf(f, "track2=fname2, \n");
  	fprintf(f, "pc=pc_fname, \n");
  	fprintf(f, "name=\"%s\", \n", fname);
  	fprintf(f, "window=\"%i\", \n", wSize);
  	fprintf(f, "kernel=\"%s\",\n", getKernelType());
  	fprintf(f, "nFgr=\"%i\",\n", nFg);
  	fprintf(f, "nBkg=\"%i\",\n", nBkg);
  	fprintf(f, "Bkg_av=\"%.4f\",\n", avBg);
  	fprintf(f, "Fg_av=\"%.4f\",\n", avFg);
  	fprintf(f, "Bkg_sd=\"%.4f\", \n", sdBg);
  	fprintf(f, "Fg_sd=\"%.4f\",\n", sdFg);
  	fprintf(f, "tot_cor=\"%.4f\",\n", totCorr);
  	fprintf(f, "Mann_Z=\"%.4f\",  \n", MannW->z);
  	fprintf(f, "p_value=\"%.2e\" \n", MannW->pVal);	
  	fprintf(f, "), output_file = file.path(getwd(), \"%s.html\"))\n", fname);
	
	fclose(f);
}

void printR(){
	char *s, b[2048], fname[1024];

	strcat(strcpy(b,outFile),".r");
	FILE *f=xopen(b,"wt");

	strcpy(b,outFile);
	s=strrchr(b,'/'); if(s==0) s=outFile; else s++; strcpy(fname,s);
//==============================================================================
	fprintf(f," #  Read the data  \n");
	fprintf(f," name <-  \'%s\'  \n\n",fname);
	fprintf(f," fg <- read.table(paste(name, '.fg', sep = '')) \n");
	fprintf(f," bkg<- read.table(paste(name, '.bkg', sep = '')) \n");
	if(writeDistCorr)
		fprintf(f," dist <- read.table(paste(name, '.dist', sep = ''), header=TRUE) \n\n");
	fprintf(f," #  Define plot limits \n\n");
	if(writeDistr==DISTR_SHORT)
		fprintf(f," y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,1])$y)) \n");
	if(writeDistr==DISTR_DETAIL)
		fprintf(f," y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y)) \n");
deb(crossWidth);
	fprintf(f," x_lim2 <- c(-%i,%i) \n\n",crossWidth,crossWidth);
	fprintf(f," # set x scale to kilobases \n");
	fprintf(f," x_lim2 <- x_lim2/1000 \n\n");
	if(writeDistr==DISTR_DETAIL){
		fprintf(f," #get chromosome data for plots, example for chr1.  \n");
		fprintf(f," #Some times you should also reset y_lim for plots \n");
		fprintf(f," #fg_chrom <- fg[fg[,1]==\"chr1\",] \n");
		fprintf(f," #dist_chrom <- dist$chr1 \n\n\n");
	}
	fprintf(f," # save plot to pdf \n");
	if(writePDF) fprintf(f," pdf(paste(name,'.pdf', sep=''), height = 9, width = 5) \n\n");
	fprintf(f," #  create the plot \n");
	fprintf(f," old.par <- par( no.readonly = TRUE ) \n");
	int nPar=1; if(writeDistCorr) nPar++;  if(doAutoCorr) nPar++; if(LCExists) nPar++;
	fprintf(f," par( mfrow = c( %i, 1 ), oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0)) \n\n",
			nPar);
	char sub[1024]; if(writePDF) sub[0]=0; else snprintf(sub,sizeof(sub),"\\n%s",fname);
	fprintf(f," plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), xlab='correlation coefficient',ylab='density', \n");
	fprintf(f," col='red', main='Distribution of correlations%s', \n",sub);
	fprintf(f," cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2) \n");
	if(writeDistr==DISTR_SHORT)
		fprintf(f," lines(density(fg[[1]]), col='blue', lwd=2) \n");
	else{
		fprintf(f," lines(density(fg[,4]), col='blue', lwd=2) \n");
	}
	fprintf(f," #plot line for chomosome \n");
	fprintf(f," #lines(density(fg[,4]), col='green', lwd=2) \n\n\n");

	if(writeDistCorr){
		fprintf(f," plot(dist$x/1000, dist$Fg, type='l',col='blue', xlim=x_lim2, \n");
		const char*dens="density";
		if(XYCorrScale!=1) sprintf(b,"%s*%i",dens,XYCorrScale); else strcpy(b,dens);
		fprintf(f," main='Cross-correlation function',xlab='Distance (kb)',ylab='%s',cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2) \n",b);
		fprintf(f," lines(dist$x/1000,dist$Bkg , col='red',lwd=2) \n");
		fprintf(f," #plot line for chomosome \n");
		fprintf(f," #lines(dist$x/1000, dist_chrom , col='green',lwd=2) \n\n");
	}
//======================================== What does it mean ???????????????????????????
//	if(LCExists){
//		const char *ss="";
//		if(LCScale==LOG_SCALE) ss="(log)";
//		fprintf(f,"\nlc=read.table(paste(name, \'.LChist\', sep = \'\'),header=TRUE)\n");
//		fprintf(f,"xmax=%.0f; xmin=%.0f; x0=xmax-150; x1=x0+20; y0=90; dy=7;\n",1000.,0.);
//		fprintf(f,"plot(lc$val, lc$r_CDF_obs*1000, col=\'blue\', lwd=2, \n");
//		fprintf(f,"type=\'l\',ylim=c(0,100), ylab=\'FDR\',");
//		fprintf(f,"xlab=\'normalized LC-values %s\',\n",ss);
//		fprintf(f,"main=\'Local correlation distributions\');\n");
//		fprintf(f,"lines(lc$val, lc$l_CDF_obs*1000, col=\'blue\', lwd=2);\n");
//		fprintf(f,"lines(lc$val, lc$l_CDF_exp*1000, col=\'red\', lwd=2);\n");
//		fprintf(f,"lines(lc$val, pmin(lc$L_FDR, lc$R_FDR), col='black', lwd=2);\n");
//
//		fprintf(f,"segments(-1000,5,1000,5,col='gray',lwd=1);");
//	}
	if(doAutoCorr){
		char bqx[4096];

		fprintf(f,"nameAuto1 <- \'%s.auto\'\n",getFnameWithoutExt(bqx, track1->name));
		fprintf(f,"nameAuto2 <- \'%s.auto\'\n",getFnameWithoutExt(bqx, track2->name));
		fprintf(f,"Auto1<- read.table(nameAuto1, sep = \'\')\n");
		fprintf(f,"Auto2<- read.table(nameAuto2, sep = \'\')\n");

		fprintf(f,"ylim1=min(min(Auto1$V2), min(Auto2$V2))\n");
		fprintf(f,"ylim2=max(max(Auto1$V2), max(Auto2$V2))\n");
		fprintf(f,"x_lim3=c(max(min(auto1$V1/1000),-100), min(max(auto1$V1/1000),100))\n");
		fprintf(f," plot  (Auto1$V1/1000,Auto1$V2,col=\'blue\',type=\'l\', xlim=x_lim3,ylim=c(ylim1,ylim2), main='Autocorrelation', xlab='distance (kb)', ylab='autocorr',lwd=2)\n");
		fprintf(f,"lines (Auto2$V1/1000,Auto2$V2,col=\'red\',lwd=2)\n");
		fprintf(f,"dy=ylim2-ylim1;  y1=ylim1+dy*0.5;  y2=ylim1+dy*0.75\n");
		fprintf(f,"text(0,y1,nameAuto1,col=\'blue\',pos = 4)\n");
		fprintf(f,"text(0,y2,nameAuto2,col=\'red\',pos = 4)\n");

	}
	fprintf(f," par( old.par ) \n\n");
	if(writePDF) fprintf(f," dev.off() \n");

	fclose(f);
}


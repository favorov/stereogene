/*
 * out.cpp
 *
 *  Created on: 12.12.2013
 *      Author: Mironov
 */
#include "track_util.h"
#include <sys/file.h>

void printR(int type);

char *r_pdf=reportPDF, *r_html=reportHTML;
statTest *MannW;
int XYCorrScale=100;
//================================================================= write results in ENCODE Broad Peak format
void printCorrelations(){
	verb("\nWrite correlations...\n");
	char b[TBS];


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
	track1->writeAuto(track1->path);
	track2->writeAuto(track1->path);
}


void printChrDistances(char *fname){
	FILE *f=xopen(fname,"wt");
	XYbgcorrelation.normilize();
	XYfgCorrelation.normilize();
	normChromDist();
	fprintf(f,"# %s vs %s\n",track1->name, track2->name);


	//========================================
	fprintf(f,"dist");
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


void printStat(){
	verb("Write statistics\n");
	writeLog("Write statistics\n");
	char b[TBS];
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
		sprintf(b,"%s.tsv",statFileName);
		fg=fileExists(b);
		if((outRes & TAB)!=0) {
			f=gopen(b,"a+t");
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
		sprintf(b,"%s.tsv",paramsFileName);
		fg=fileExists(b);
		if((outRes&TAB)!=0){
			f=gopen(b,"a+t");
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
			snprintf(b,sizeof(b), "%s.xml",statFileName);
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
		fprintf(f,"%9.2f\t%9.3f\t%9.3f\n",l,spectrumX[i],spectrumY[i]);
	}


	fclose(f);
}


//============================================== write the foreground and background distributions of the correlations
void printBgDistr(){
	char b[TBS];
	strcat(strcpy(b,outFile),".bkg");					// open file for background observations
	FILE* fbkg=xopen(b,"wt");
	for(int i=0; i<nBkg; i++) fprintf(fbkg,"%f\n",BkgSet[i]);
	fclose(fbkg);
}


void printFgDistr(){
	char b[TBS];
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
//const char *template1="report_r_template1.Rmd";
//const char *template2="report_r_template2.Rmd";
//const char *template3="report_r_template3.Rmd";







void printRaw(FILE *f, const char * prm, double val){
	if(val >0.1)
		fprintf(f,"	<TR VALIGN=TOP>	<TD> %s </TD> <TD>	%.2f	</TD>	</TR>\n",prm,val);
	else
		fprintf(f,"	<TR VALIGN=TOP>	<TD> %s </TD> <TD>	%.2e	</TD>	</TR>\n",prm,val);
}
void printRaw(FILE *f, const char * prm, int val){
	fprintf(f,"	<TR VALIGN=TOP>	<TD > %s </TD> <TD>	%i	</TD>	</TR>\n",prm,val);
}
void printRaw(FILE *f, const char * prm, const char * val){
	fprintf(f,"	<TR VALIGN=TOP>	<TD > %s </TD> <TD>	%s	</TD>	</TR>\n",prm,val);
}

const char*pdftext="  text(x=0,y=y ,adj = c(0,1), cex=0.8, label='";
void printPDFRaw(FILE *f, const char * prm, double val){
	if(val >0.1)
		fprintf(f,"%s %s = %.2f'); y=y-dy\n",pdftext,prm,val);
	else
		fprintf(f,"%s %s = %.2e'); y=y-dy\n",pdftext,prm,val);
}
void printPDFRaw(FILE *f, const char * prm, int val){
	fprintf(f,"%s %s = %i'); y=y-dy\n",pdftext,prm,val);
}
void printPDFRaw(FILE *f, const char * prm, const char * val){
	fprintf(f,"%s %s = %s'); y=y-dy\n",pdftext,prm,val);
}



void printHTML(){
	char fname[TBS];

	snprintf(fname,sizeof(fname),"%s%s.html",curRepPath,curOutFname);

	FILE *f=xopen(fname,"w");


	fprintf(f,"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0\">\n");
	fprintf(f,"<HTML>\n<BODY>\n");


	fprintf(f,"<TABLE>\n	<TR VALIGN=TOP>\n"); //begin tab tables


	fprintf(f,"<TD>\n<H2> Input </H2>\n<TABLE BORDER=1 BORDERCOLOR=\"#000000\" CELLPADDING=4 CELLSPACING=0>\n");//begin params table
	printRaw(f,"Track1",track1->name);
	printRaw(f,"Track2",track2->name);
	printRaw(f,"Window Size (kb)",wSize/1000);
	printRaw(f,"Bin Size"      ,binSize);
	printRaw(f,"Kernel width"  ,kernelSigma);
	printRaw(f,"Max zero (%)"  ,maxZero0);
	printRaw(f,"Max NA (%)"    ,maxNA0);
	fprintf(f,"</TABLE>\n");//end param table


	fprintf(f,"<TD>\n<H2> Results </H2>\n<TABLE BORDER=1 BORDERCOLOR=\"#000000\" CELLPADDING=4 CELLSPACING=0>\n");//begin res table
	printRaw(f,"Output files",outFile);
	printRaw(f,"n Foreground comparisons"	,nFg);
	printRaw(f,"n Background comparisons"	,nBkg);
//	printRaw(f,"Foregr. total correlation"  ,totCorr);
	printRaw(f,"Foregr. correlation",avFg);
	printRaw(f,"Foregr. correlation std.dev",sdFg);
//	printRaw(f,"Backgr. total correlation"  ,BgTotal);
	printRaw(f,"Backgr. correlation",avBg);
	printRaw(f,"Backgr. correlation std.dev",sdBg);
	printRaw(f,"Mann-Witney p-value"    ,mannW_p);
	fprintf(f,"</TABLE>\n");//end param table


	fprintf(f,"</TABLE><hr/>\n");//end tab tables

	//================================= Make svg file
	fprintf(f,"<img src=\"%s.svg\" alt=\"???\" width=70%% ALIGN=\"left\">\n", curOutFname);
	fprintf(f,"</BODY>\n</HTML>");
}


void printR(){
	reportPDF [0]=0;
	reportHTML[0]=0;

	if(RScriptFg & R) 	{		printR(R);		}
	if(RScriptFg & PDF) {		printR(PDF);					strcat(strcpy(reportPDF ,curOutFname), ".pdf" );}
	if(RScriptFg & HTML){		printR(HTML);	printHTML();	strcat(strcpy(reportHTML,curOutFname), ".html");}
}


void printR(int type){

	char b[4096], fname[4096],rFile[4096];
	int pH=plotH;
	if(type==R)		strcat(strcpy(rFile,outFile),".r");
	if(type==PDF)	strcat(strcpy(rFile,outFile),"_pdf.r");
	if(type==HTML)	strcat(strcpy(rFile,outFile),"_svg.r");
	FILE *f=xopen(rFile,"wt");

	getFnameWithoutPath(fname, outFile);


//==============================================================================
	fprintf(f," #  Read the data  \n");
	fprintf(f," name <-  \'%s\'  \n\n",fname);
	if(reportPath){
		fprintf(f," report <-  \'%s%s\'  \n\n",reportPath,fname);
	}
	else{
		fprintf(f," report <-  \'%s\'  \n\n",fname);
	}
	fprintf(f," fg <- read.table(paste(name, '.fg', sep = '')) \n");
	fprintf(f," bkg<- read.table(paste(name, '.bkg', sep = '')) \n");
	if(writeDistCorr){
		fprintf(f," dist <- read.table(paste(name, '.dist', sep = ''), header=TRUE) \n\n");
	}

	fprintf(f," #  Define plot limits \n\n");
	if(writeDistr==DISTR_SHORT){
		fprintf(f," y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,1])$y)) \n");
	}
	if(doAutoCorr || writeDistr==DISTR_DETAIL) pH+=plotH;
	if(writeDistr==DISTR_DETAIL)
		fprintf(f," y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y)) \n");
	fprintf(f," x_lim2 <- c(-%i,%i) \n\n",crossWidth,crossWidth);
	fprintf(f," # set x scale to kilobases \n");
	fprintf(f," x_lim2 <- x_lim2/1000 \n\n");
	if(writeDistr==DISTR_DETAIL){
		fprintf(f," #get chromosome data for plots, example for chr1.  \n");
		fprintf(f," #Some times you should also reset y_lim for plots \n");
		fprintf(f," #fg_chrom <- fg[fg[,1]==\"chr1\",] \n");
		fprintf(f," #dist_chrom <- dist$chr1 \n\n\n");
	}
	fprintf(f," # save plot to pdf/SVG \n");

	if(type== PDF)  fprintf(f,"pdf(paste(report,'.pdf', sep=''), height = %i, width = %i) \n\n",pH+plotH,plotW);
	if(type== HTML) fprintf(f,"svg(paste(report,'.svg', sep=''), height = %i, width = %i) \n\n",pH,plotW);
	fprintf(f," #  create the plot \n");
	fprintf(f," old.par <- par( no.readonly = TRUE ) \n");
	int nPar=1;
	if(doAutoCorr || writeDistr==DISTR_DETAIL) nPar++;
	if(type== PDF)  nPar++;
	fprintf(f," par( mfrow = c( %i, 2 ), oma = c( 0, 0, 0, 0 ),mar=c(3,3,3,1),mgp=c(1.6,0.45,0)) \n\n",
			nPar);

	//	===============================================================================================
	//	=============================== print pdf text ================================================
	//	===============================================================================================
	if(type==PDF){
	fprintf(f,"plot(x = 0:1, y = 0:1,bty = 'n',type = 'n', xlab='',ylab='',\n");
	fprintf(f,"      xaxt = 'n', yaxt = 'n', xlim=c(0,1), ylim=c(0,2), main='Input')\n\n");
	fprintf(f," y=2; dy=0.2\n\n");
	printPDFRaw(f,"Track1",track1->name);
	printPDFRaw(f,"Track2",track2->name);
	printPDFRaw(f,"Window Size (kb)",wSize/1000);
	printPDFRaw(f,"Kernel width"    ,kernelSigma);
	printPDFRaw(f,"Bin Size"        ,binSize);
	printPDFRaw(f,"Max zero (%)"    ,maxZero0);
	printPDFRaw(f,"Max NA (%)"      ,maxNA0);

	//	===============================================================================================
	fprintf(f,"plot(x = 0:1, y = 0:1,bty = 'n',type = 'n', xlab='',ylab='',\n");
	fprintf(f,"      xaxt = 'n', yaxt = 'n', xlim=c(0,1), ylim=c(0,2), main='Results')\n\n");
	fprintf(f," y=2; dy=0.2\n\n");
	printPDFRaw(f,"Output files",fname);
	printPDFRaw(f,"n Foreground comparisons",nFg);
	printPDFRaw(f,"n Background comparisons",nBkg);
//	printPDFRaw(f,"Foregr. total correlation"       ,totCorr);
	printPDFRaw(f,"Foregr.  correlation"     ,avFg);
	printPDFRaw(f,"Foregr. correlation std.dev"   ,sdFg);
//	printPDFRaw(f,"Backgr. total correlation"       ,BgTotal);
	printPDFRaw(f,"Backgr.  correlation"     ,avBg);
	printPDFRaw(f,"Backgr. correlation std.dev"   ,sdBg);
	printPDFRaw(f,"Mann-Witney p-value"    ,mannW_p);
	}

	//	===============================================================================================
	//	===============================================================================================
	//	===============================================================================================
	char sub[1024];
	snprintf(sub,sizeof(sub),"\\n%s",fname);

	const char* cex="      cex.axis = 0.8,  cex.lab = 1,  cex.main = 1,lwd=2) \n";
	fprintf(f," y_lim2 <- max(max(dist$Fg),max(dist$Bkg))\n");
	fprintf(f," plot(density(bkg[[1]]),main='Distribution of correlations%s',\n",sub);
			fprintf(f,	"      xlim=c(-1,1), ylim=c(0, y_lim1),\n");
			fprintf(f,	"      xlab='correlation coefficient',ylab='density',  col='red', \n");
			fprintf(f,"%s",cex);
	fprintf(f," legend(-1, y_lim1, legend=c('Foreground','Background'),\n");
    fprintf(f,"     col=c('blue','red'), lty=1:2, cex=0.5)\n");
	if(writeDistr==DISTR_SHORT)
		fprintf(f," lines(density(fg[[1]]), col='blue', lwd=2) \n");
	else{
		fprintf(f," lines(density(fg[,4]), col='blue', lwd=2) \n");
	}
	fprintf(f," #plot line for chomosome \n");
	fprintf(f," #lines(density(fg[,4]), col='green', lwd=2) \n\n\n");


	const char*dens="density";
	if(writeDistCorr){
		if(XYCorrScale!=1) snprintf(b,sizeof(b), "%s*%i",dens,XYCorrScale); else strcpy(b,dens);
		fprintf(f," plot(dist$dist/1000, dist$Fg, type='l', main='Cross-correlation function%s',\n",sub);
		fprintf(f,"      xlim=x_lim2, xlab='Distance (kb)',ylab='%s',col='blue',\n",b);
		fprintf(f,"%s",cex);
		fprintf(f," legend(x_lim2[1], y_lim2, legend=c('Foreground','Background'),\n");
		fprintf(f,"     col=c('blue','red'), lty=1:2, cex=0.5)\n");

		fprintf(f," lines(dist$dist/1000,dist$Bkg , col='red',lwd=2) \n");
		fprintf(f," #plot line for chomosome \n");
		fprintf(f," #lines(dist$x/1000, dist_chrom , col='green',lwd=2) \n\n");
	}
//============================================ Plot for chromosomes
	if(writeDistr==DISTR_DETAIL){
		fprintf(f,"#================= Plot distrib. by chromosomes\n");
		fprintf(f,"chrDist=list(");
		for(int i=0; i<n_chrom; i++){
			if(i) fprintf(f,",");
			fprintf(f,"       fg[fg$V1==\'%s\',]$V4\n", chrom_list[i].chrom);
		}
		fprintf(f,")\n");
		fprintf(f,"names=c(");
		for(int i=0; i<n_chrom; i++){
			if(i) fprintf(f,",");
			fprintf(f,"       \'%s'\n", chrom_list[i].chrom);
		}
		fprintf(f,")\n");
		fprintf(f,"boxplot(chrDist,names=names,las = 2,cex.lab = 1,cex.axis=0.8)\n");

	}

	if(doAutoCorr){
		fprintf(f,"#================= Plot autocorr\n");
		fprintf(f,"nameAuto1 <- \'%s.auto\'\n", track1->name);
		fprintf(f,"nameAuto2 <- \'%s.auto\'\n", track2->name);
		fprintf(f,"Auto1<- read.table(nameAuto1, sep = \'\')\n");
		fprintf(f,"Auto2<- read.table(nameAuto2, sep = \'\')\n");
		fprintf(f,"y_lim3 <- max(max(Auto1$V2),max(Auto2$V2))\n");


		fprintf(f,"ylim1=min(min(Auto1$V2), min(Auto2$V2))\n");
		fprintf(f,"ylim2=max(max(Auto1$V2), max(Auto2$V2))\n");
		fprintf(f," plot  (Auto1$V1/1000,Auto1$V2,type=\'l\', main='Autocorrelation\\n%s,\\n%s',\n",track1->name,track2->name);
		fprintf(f,"      xlim=x_lim2, ylim=c(ylim1,ylim2),  xlab='distance (kb)', ylab='autocorr',col=\'blue\',\n");
		fprintf(f,"%s",cex);
		fprintf(f," lines  (Auto2$V1/1000,Auto2$V2,\n");
		fprintf(f,"      xlim=x_lim2, ylim=c(ylim1,ylim2),  xlab='distance (kb)', ylab='autocorr',col=\'red\',\n");
		fprintf(f,"%s",cex);
		fprintf(f," legend(x_lim2[1], 1, legend=c('%s','%s'),\n",track1->name,track2->name);
	    fprintf(f,"     col=c('blue','red'), lty=1:2, cex=0.5)\n");


	}
	fprintf(f," par( old.par ) \n\n");
	if(type==PDF || type==	HTML) fprintf(f,"%s"," dev.off() \n");

	fclose(f);

	if(Rscript && (type==PDF || type==	HTML)){
		char b[4096], cwd[4096];
		char *rf=strrchr(rFile,'/');
		if(rf) rf=rf+1;
		else   rf=rFile;
		getcwd(cwd, sizeof(cwd));
		strcpy(b,outFile);
		char *s=strrchr(b,'/');
		if(s) *s=0;
		chdir(b);
		snprintf(b,sizeof(b), "\"%s\" %s --vanilla",Rscript,rf);
		system(b);
		chdir(cwd);
	}
}

//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================

void printParamNames(FILE* f){
	fprintf(f,"id\tversion");
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg && (pparams[i]->prog&progType)) fprintf(f,"\t%s",pparams[i]->name);
	}
	fprintf(f,"\n");
}
NamedRes *results[]={
		new NamedRes("id",printId),
		new NamedRes("Date",dateTime),
		new NamedRes("version", (char**)&version),

		new NamedRes("input"),
		new NamedRes("name1",&trackName1),
		new NamedRes("name2",&trackName2),

		new NamedRes("res"),
		new NamedRes("nFgr",&nFg),
		new NamedRes("nBkg",&nBkg),
		new NamedRes("Fg. Corr",&avFg),
		new NamedRes("FgCorr_sd",&sdFg),
		new NamedRes("Bg. Corr",&avBg),
		new NamedRes("BgCorr_sd",&sdBg),
		new NamedRes("Mann-Z",&mannW_Z),
		new NamedRes("p-value",&mannW_p),
		new NamedRes("PDF_report" ,&r_pdf),
		new NamedRes("HTML_report",&r_html),
		0
};


void printStat(FILE *f){
	for(int i=0; results[i]; i++){
		if(results[i]->type==0) continue;
		if(i)fprintf(f,"\t");
		results[i]->printValue(f);
	}
	for(int i=0; pparams[i] ; i++){
		if((pparams[i]->printFg&2)==2){
			if(i)fprintf(f,"\t");
			pparams[i]->printParamValue(f);
		}
	}
	fprintf(f,"\n");
}


void printParams(FILE* f){
	char b[4096];
	fprintf(f,"%08lx\t%s",id,version);
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg && (pparams[i]->prog&progType))
			fprintf(f,"\t%s",pparams[i]->printParamValue(b));
	}
	fprintf(f,"\n");
}
void printXMLparams(FILE *f){
	char b[4096];
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg && (pparams[i]->prog&progType))
			fprintf(f,"%s=\"%s\" ",pparams[i]->name,pparams[i]->printParamValue(b));
	}
}
void printStatHeader(FILE *f){
	for(int i=0; results[i]; i++){
		if(results[i]->type==0) continue;
		if(i)fprintf(f,"\t");
		fprintf(f,"%s",results[i]->name);
	}
	for(int i=0; pparams[i] ; i++){
		if((pparams[i]->printFg&2)==2)
			fprintf(f,"\t%s",pparams[i]->name);
	}
	fprintf(f,"\n");
}



void printXML(FILE *f){
	char b[1024];
	fprintf(f,"<run ");
	int kk=0;
	for(int i=0; results[i]; i++){
		if(results[i]->type==0) {
			if(kk) fprintf(f,"/>\n\t");
			else fprintf(f,">\n\t");
			kk++;
			fprintf(f,"<%s ",results[i]->name); continue;
		}
		fprintf(f," %s=\"%s\"",results[i]->name,results[i]->printValue(b));
	}
	fprintf(f,"/>\n\t<prm ");
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg && (pparams[i]->prog&progType))
			fprintf(f,"%s=\"%s\" ",pparams[i]->name, pparams[i]->printParamValue(b));
	}
	fprintf(f,"/>\n</run>\n");
}






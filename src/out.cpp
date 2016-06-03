/*
 * out.cpp
 *
 *  Created on: 12.12.2013
 *      Author: 1
 */
#include "track_util.h"
#include <sys/file.h>

statTest *MannW;
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
	if(outSpectr) correlation.printSpect(strcat(strcpy(b,outFile),".spect"));
	if(outChrom ) printChomosomes(strcat(strcpy(b,outFile),".chrom"));
	errStatus=0;
}

void printChrDistances(char *fname){
	FILE *f=xopen(fname,"wt");
	bgcorrelation.norm();
	correlation.norm();
	fprintf(f,"x\tBkg\tFg\tFgPlus\tFgMinus");

	if(writeDistCorr==CHR_DETAIL){
		for(int i=0; i<n_chrom; i++){
			Chromosome *chr=chrom_list+i;
			if(chr->densCount)	fprintf(f,"\t%s",chr->chrom);
		}
	}
	fprintf(f,"\n");
	for(int j=0; j<profWithFlanksLength; j++){
		int k=(j+profWithFlanksLength/2)%profWithFlanksLength;
		fprintf(f,"%i",(j-profWithFlanksLength/2)*stepSize);
		fprintf(f,"\t%9.5f\t%9.5f\t%9.5f\t%9.5f", bgcorrelation.correlation[k]*100,
				correlation.correlation[k]*100, correlation.corrPlus[k]*100, correlation.corrMinus[k]*100);
		if(writeDistCorr==CHR_DETAIL){
			for(int i=0; i<n_chrom; i++){
				Chromosome *chr=chrom_list+i;
				if(chr->densCount){
					if(chr->count > 1){
						double e=bTrack1.av*bTrack2.av;
						double d=(bTrack1.sd*bTrack2.sd*profWithFlanksLength*(chr->count));    //*chr->count
						double x=(chr->distDens[k]-e*chr->count)/d;

						fprintf(f,"\t%7.3f",x*100);
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

void printStat(){
	verb("Write statistics\n");
	char b[2048];

	MannW=MannWhitney(FgSet, nFg, BkgSet, nBkg);	// do Mann-Whitney test
	xverb("p-val=%e nWindows=%i Kern=%s\n",
		MannW->pVal, nFg, getKernelType());

	double avBg,sdBg,avFg,sdFg;
	getStat(FgSet,nFg,avFg,sdFg);
	getStat(BkgSet,nBkg,avBg,sdBg);

	FILE *f=0;
	const char *pcname = "-";
	char pc = '-';
	if (pcorProfile!=0){
		pc = '+';
		pcname = pcorProfile;
	}

	bool fg=fileExists(statFileName);
	if((outRes & TAB)!=0) {
		f=gopen(statFileName,"a+t");
		if(f!=0) flockFile(f);
		else {
			fprintf(stderr,"Can not open file %s\n",statFileName);
			writeLog("Can not open file %s\n",statFileName);
		}
	}
	if(!fg && f){								//================ write the header
		fprintf(f,"%-6s\t%-20s\t%-20s","id","name1","name2");
		fprintf(f,"\t%-6s\t%-6s\t%-6s\t%-6s","window","kern","nFgr","nBkg");
		fprintf(f,"\t%-8s\t%-8s\t%-8s\t\%-8s","Bkg_av","Fg_av","Bkg_sd","Fg_sd");
		fprintf(f,"\t%-9s\t%-8s\t%-8s\t%-7s","tot_cor","avCorr", "Mann-Z","p-value");
		fprintf(f,"\t%-6s\n", "pc");
	}
	//==================================================== write the statistics
	if(f){
		fprintf(f,"%lx\t%-10s\t%-10s",id,alTable.convert(bTrack1.name), alTable.convert(bTrack2.name));
		fprintf(f,"\t%-6i\t%-6s\t%6i\t%6i",wSize,getKernelType(),nFg, nBkg);
		fprintf(f,"\t%8.4f\t%8.4f\t%8.4f\t%8.4f",avBg,avFg,sdBg,sdFg);
		fprintf(f,"\t%8.4f\t%8.4f\t%8.1f\t%-7.2e\t%-6c\n",totCorr,FgAvCorr,MannW->z,MannW->pVal, pc);
		funlockFile(f);
		fclose(f);
	}

	//================================================== write parameters
	writeLog("Write params...\n");
	fg=fileExists(paramsFileName);
	if((outRes&TAB)!=0){
		f=gopen(paramsFileName,"a+t");
		if(f) flockFile(f);
		else {
			fprintf(stderr,"Can not open file %s\n",paramsFileName);
			writeLog("Can not open file %s\n",paramsFileName);
		}
	}
	if(!fg && f){								//================ write the header
		fprintf(f,"%-6s\t%-20s\t%-20s","id","trackPath","resPath");
		fprintf(f,"\t%-20s\t%-12s\t%-20s","map","mapIv","pcorProfile");
		fprintf(f,"\t%-2s\t%-6s\t%-6s","NA", "maxNA","maxZer");
		fprintf(f,"\t%-6s\t%-6s\t%-6s","interv", "strand","compl");
		fprintf(f,"\t%-4s\t%-6s","step", "bpType");
		fprintf(f,"\t%-6s\t%-6s\t%-6s\t%-6s","wSize", "wStep","flank","noise");
		fprintf(f,"\t%-6s\t%-8s\t%-8s","kernel", "Kern-Sgm","kern-Sh");
		fprintf(f,"\t%-8s\t%-8s\t%-8s\n","nShuffle", "MaxShfl","threshold");
	}
	//==================================================== write the parameters
	if(f){
		fprintf(f,"%lx\t%-20s\t%-20s",id,trackPath,resPath);
		const char *mf=mapFil ? mapFil : "-";
		fprintf(f,"\t%-20s\t%-12s\t%-20s",mf,miv.print(b),pcname);
		fprintf(f,"\t%-2i\t-%6.1f\t%-6.1f",NAFlag, maxNA0,maxZero0);
		fprintf(f,"\t%-6i\t%-6i\t%-6i",intervFlag0, strandFg0,complFg);
		fprintf(f,"\t%-4i\t%-6i",stepSize, bpType);
		fprintf(f,"\t%-6i\t%-6i\t%-6i\t%-6.2f",wSize, wStep,flankSize,noiseLevel);
		fprintf(f,"\t%-6s\t%-8.0f\t%-8.0f",getKernelType(), kernelSigma,kernelShift);
		fprintf(f,"\t%-i\t%-i\t%-8i\n",nShuffle, maxShuffle,threshold);
		funlockFile(f);
		fclose(f);
	}
	writeLog("OK\n");

	writeLog("Write XML ...\n");
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
		fprintf(xml,"<run id=\"%lx\" ver=\"%s\">\n", id, version);
		fprintf(xml,"\t<input track1=\"%s\" track2=\"%s\"/>\n",bTrack1.name,bTrack2.name);
		fprintf(xml,"\t<output out=\"%s\"/>\n",outFile);
		fprintf(xml,"\t<prm ");
		fprintf(xml,"kernel=\"%s\" ",getKernelType());
		fprintf(xml,"kernelSigma=\"%.0f\" ",kernelSigma);
		fprintf(xml,"kernelShift=\"%.0f\" ",kernelShift);
		fprintf(xml,"kernelNS=\"%.0f\" ",kernelNS*100);
		fprintf(xml,"wSize=\"%i\" ",wSize);
		fprintf(xml,"wStep=\"%i\" ", wStep);
		fprintf(xml,"flankSize=\"%i\" ",flankSize);
		fprintf(xml,"noiseLevel=\"%.2f\" ",noiseLevel);
		fprintf(xml,"nShuffle=\"%i\" ",nShuffle);
		fprintf(xml,"maxShuffle=\"%i\" ", maxShuffle);
		fprintf(xml,"threshold=\"%i\" ",threshold);
		fprintf(xml,"stepSize=\"%i\" ",stepSize);
		fprintf(xml,"bpType=\"%i\" ", bpType);
		fprintf(xml,"NAFlag=\"%i\" ",NAFlag);
		fprintf(xml,"maxNA=\"%.1f\" ", maxNA0);
		fprintf(xml,"maxZero=\"%.1f\" ",maxZero0);
		fprintf(xml,"intervFlag=\"%x\" ",intervFlag0);
		fprintf(xml,"strandFg=\"%i\" ", strandFg0);
		fprintf(xml,"complFg=\"%i\" ",complFg);
		if(mapFil) fprintf(xml,"mapFile=\"%s\" miv=\"%s\"",mapFil,b);
		if (pcorProfile!=0) fprintf(xml,"pcname=\"%s\" ",pcname);
		fprintf(xml,"/>\n");

		fprintf(xml,"\t<res ");
		fprintf(xml,"nFg=\"%i\" ",nFg);
		fprintf(xml,"nBkg=\"%i\" ", nBkg);
		fprintf(xml,"totCorr=\"%.4f\" ", totCorr);
		fprintf(xml,"FgAvCorr=\"%.4f\" ", FgAvCorr);
		fprintf(xml,"BgAvCorr=\"%.4f\" ", BgAvCorr);
		fprintf(xml,"MannZ=\"%.4f\" ",MannW->z);
		fprintf(xml,"pVal=\"%.2e\" ",MannW->pVal);

		fprintf(xml,"avBg=\"%.4f\" ",avBg);
		fprintf(xml,"avFg=\"%.4f\" ",avFg);
		fprintf(xml,"sdBg=\"%.4f\" ",sdBg);
		fprintf(xml,"sdFg=\"%.4f\" ",sdFg);
		fprintf(xml,"/>\n");
		fprintf(xml,"</run>\n");
		fprintf(xml,"</xml>\n");
		funlockFile(xml);
		fclose(xml);
	}
	writeLog("WriteStat - OK\n");
}

void Correlation::printSpect(char *fname){
	FILE *f=xopen(fname,"wt");
	double dx=0,dy=0;
	for(int i=0; i<profWithFlanksLength; i++){
		dx+=spectrumX[i]; dy+=spectrumY[i];
	}
	for(int i=0; i<profWithFlanksLength; i++){
		spectrumX[i]=sqrt(spectrumX[i]/dx*profWithFlanksLength);
		spectrumY[i]=sqrt(spectrumY[i]/dy*profWithFlanksLength);
	}

	for(int i=5; i<profWithFlanksLength/2; i++){
		double l=(double)(profWithFlanksLength)*stepSize/(i+1);
		fprintf(f,"%.2f\t%g\t%g\n",l,spectrumX[i],spectrumY[i]);
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
	fclose(fBpeak);
	xfree(pPair);
}



void printR(){
	writeLog("Write R ...");
	char *fn=alTable.convert(outFile), *s,b[2048], fname[1024];
	const char *cex="cex.axis = 0.8,  cex.lab = 0.8,  cex.main = 0.8", *lwd="lwd=2",
			*lab="xlab=\'correlation coefficient\',ylab=\'density\'";

	int    x0,x1;
	double y0,y1;
	correlation.getLimits(x0,x1, y0, y1);

	strcat(strcpy(b,fn),".r");
	FILE *f=xopen(b,"wt");

	strcpy(b,outFile);
	s=strrchr(b,'/'); if(s==0) s=outFile; else s++; strcpy(fname,s);

	fprintf(f,"#  Read the data \n\n");

	fprintf(f,"fg <- read.table(\'%s.fg\')\n",fname);
	fprintf(f,"bkg<- read.table(\'%s.bkg\')\n",fname);

	fprintf(f,"dist <- read.table(\'%s.dist\', header=TRUE)\n",fname);
	fprintf(f,"\n#  Define plot limits\n\n");

	fprintf(f,"y_lim1 <- max(max(density(bkg[,1])$y),max(density(fg[,4])$y))\n");
	fprintf(f,"y_lim2 <- c(%.2f,%.2f)\n",y0,y1);

	fprintf(f,"x_lim2 <- c(%i,%i)\n\n",x0,x1);

	fprintf(f,"# set x scale to kilobases\n");
	fprintf(f,"x_lim2 <- x_lim2/1000\n");

	fprintf(f,"\n#  create the plot\n\n");

	fprintf(f,"old.par <- par( no.readonly = TRUE )\n");
	fprintf(f,"par( mfrow = c( 2, 1 ), oma = c( 0.5, 0, 2, 0 ),mar=c(2.5,3,1.5,1),mgp=c(1,0.3,0))\n\n");

	fprintf(f,"plot(density(bkg[[1]]), xlim=c(-1,1), ylim=c(0, y_lim1), %s,\n",lab);
	fprintf(f,"col=\'red\', main=\'Distribution of correlations, p-value=%.1e\',\n",MannW->pVal );
	fprintf(f,"%s,%s)\n\n",cex,lwd);
	fprintf(f,"lines(density(fg[,4]), col=\'blue\', %s)\n\n",lwd);

	fprintf(f,"plot(dist$x/1000, dist$Fg, type=\'l\',col=\'blue\', ylim=y_lim2, xlim=x_lim2,\n");
	fprintf(f,"main=\'Cross-correlation function\',xlab=\'Distance (kb)\',ylab=\'density,%%\',%s,%s)\n",cex,lwd);
	fprintf(f,"#lines(dist$x/1000,dist$FgPlus, col=\'cyan\',%s)\n",lwd);
	fprintf(f,"#lines(dist$x/1000,dist$FgMinus, col=\'brown\',%s)\n",lwd);
	fprintf(f,"lines(dist$x/1000,dist$Bkg , col=\'red\',%s)\n",lwd);

	char fn1[1024], fn2[1024];
	strcpy(fn1,fname);
	s=strchr(fn1,'~');
	if(s){*s=0; s++; strcpy(fn2,s);}
	else {*fn2=0;}

	fprintf(f,"title(\'%s\\n%s\',cex.main = 0.8, outer = TRUE)\n\n", fn1,fn2);
	fprintf(f,"par( old.par )\n");

	fclose(f);
	writeLog("OK\n");
}


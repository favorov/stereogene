/*
 * outwig.cpp
 * Module does output of the local correlations into wig file.
 *
 *  Created on: Dec 11, 2013
 *      Author: mironov
 */
#include "track_util.h"

//TODO linear instead  log transform;
//TODO Normalize histogram;
//TODO CAlulate and print FDR


struct ScaledAray{
	float *array;		// values
	char  *counts;		// counts of observations in given points
	int curPos;			// the last position of filled part of the array
	ScaledAray();
	void init();		// initiation : allocate array. Size = (genome size) / stepSize
	void addArray(double *f, int pos);		// add array of doubles to the scaled array
	void write();							// write output files
};


ScaledAray wigCorr;
double *smoothProf2, *smoothProf1;
double minLC,maxLC;		// min and max values of the local correlation


struct WigHist:DinHistogram{
	double *FDR;
	WigHist(int l);
	void print(FILE*f);
	void fin();
};

double normWig(double v){
	return 1000.*(v-minLC)/(maxLC-minLC);
}

WigHist dHist=WigHist(500);
WigHist::WigHist(int l):DinHistogram(l){
	getMem(FDR,l,"WigHist:getFDR");
}

void WigHist::fin(){
	min =normWig( min); max= normWig( max);
	hMin=normWig(hMin); hMax=normWig(hMax);
	bin=(hMax-hMin)/l;
	DinHistogram::fin();					// normalize results
	double obs=0, exp=0;
	for(int i=l-1; i>=0; i--){
		obs+=hist[0][i]*bin;   				// observed
		exp+=hist[1][i]*bin;				// expected
		double psi=1/n[0]*bin;				// psudocounts
		double x=(exp+psi)/(obs+psi);
		FDR[i]=x;
	}
}


void WigHist::print(FILE* f){						// print the histogram
	fprintf(f,"#  min=%.3f max=%.3f \n",min,max);
	fprintf(f,"#  hMin=%.2ef  hMax=%.2ef bin=%.3f\n",hMin,hMax,bin);
	fprintf(f,"#  Observed: e0=%.3f sd0=%.3f n0=%i\n",e[0],sd[0],n[0]);
	fprintf(f,"#  Expected: e1=%.3f sd1=%.3f n1=%i\n",e[1],sd[1],n[1]);
	fprintf(f,"value\tobs\tnObs\texp\tnExp\tFDR(%%)\n");
	int i0=0, i1=l;
	for(int i=0; i<l; i++){
		if(i1==l && (hist[0][i]!=0 || hist[1][i]!=0)) i0=i;
		if(hist[0][i]!=0 || hist[1][i]!=0) i1=i;
	}

	for(int i=i0; i<=i1; i++){
		double h0=hist[0][i];
		double h1=hist[1][i];
		double fdr=FDR[i]*100;
		fprintf(f,"%.1f\t%.2e\t%6i\t%.2e\t%6i\t%.2f\n", getValue(i),h0,cnts[0][i],h1,cnts[1][i],fdr);
	}
}

//==========================================================================
ScaledAray::ScaledAray(){
	minLC=1.e+128; maxLC=-1.e+128; curPos=0; array=0; counts=0;
}

void ScaledAray::init(){
	getMem0(array,profileLength+10, "ScaledAray::init #1");
	getMem0(counts,profileLength+10, "ScaledAray::init #2");
	zeroMem(array,profileLength+10);
	zeroMem(counts,profileLength+10);
	minLC=1.e+128; maxLC=-1.e+128; curPos=0;

}


void ScaledAray::addArray(double *f, int pos){
	for(int i=0, j=pos; i<wProfSize; i++,j++){
		double x=f[i],xx=x;
		if(j<curPos){
			xx=(array[j]*counts[j]+x)/(counts[j]+1);
		}
		array[j]=xx;
		counts[j]++;
	}
	curPos=pos+wProfSize;
}



void ScaledAray::write(){
	char bf[1024];
	sprintf(bf,"%s.wig",outFile);
	FILE *outWigFile=xopen(bf,"wt");
	verb("\nwrite Local Correlations\n");
	ScoredRange pos, pos0;

	minLC=1.e+128; maxLC=-1.e+128;
	for(int i=0; i<profileLength; i++){
		if(minLC > array[i]) minLC = array[i];
		if(maxLC < array[i]) maxLC = array[i];
	}

	fprintf(outWigFile,"track type=wiggle_0 ");
	char *s=strrchr(outFile,'/'); if(s==0) s=outFile; else s++;
	fprintf(outWigFile,"description=\"correlation: %s\" \n",s);
	fprintf(outWigFile,"#Param ");
	if((outWIG & WIG_BASE)  == WIG_BASE  ) fprintf(outWigFile,"BASE "  );
	if((outWIG & WIG_CENTER)== WIG_CENTER) fprintf(outWigFile,"CENTER ");
	if((outWIG & WIG_SUM)   == WIG_SUM   ) fprintf(outWigFile,"SUM "   );
	if((outWIG & WIG_MULT)  == WIG_MULT  ) fprintf(outWigFile,"MULT "  );
	fprintf(outWigFile,"\n");
	fprintf(outWigFile,"#Scale data: min=%.4f;  max=%.4f;   scale: x=1000*(LC-min)/(max-min)\n", minLC,maxLC);
	fprintf(outWigFile,"#Source statstics: av1=%.4f;  av2=%.4f; sd1=%.4f; sd2=%.4f\n", bTrack1.av0,bTrack2.av0, bTrack1.sd0,bTrack2.sd0);

	for(int i=0; i<profileLength; i++){
		if(i%1000000 ==0) verb("%5.1f%%\r",1.*i/profileLength*100);
		int k=(int)normWig(array[i]);
		if(k>=outThreshold){
			filePos2Pos(i,&pos,binSize);
			if(pos.chrom != pos0.chrom  || pos.beg-pos0.beg > binSize){
				fprintf(outWigFile,"fixedStep chrom=%s start=%li step=%i span=%i\n",
						pos.chrom,pos.beg,binSize,binSize);
			}
			pos0.chrom=pos.chrom; pos0.beg=pos.beg;
			fprintf(outWigFile,"%i\n",k);
		}
	}
	fclose(outWigFile);
	verb("\n");
	//========================================== Write histograms ===============
	sprintf(bf,"%s.LChist",outFile);
	outWigFile=fopen(bf,"wt");
	dHist.fin();
	dHist.print(outWigFile);
	fclose(outWigFile);
	dHist.clear();
}
//==================================================================
void initOutWig(){
	if(outWIG==NONE) return;
	wigCorr.init();
}

void finOutWig(){
	if(outWIG==NONE) return;
	//====================== write correlation
	wigCorr.write();

}

//===================================================================
//============= calculate smoothed profile  c=\int f*\rho
void calcSmoothProfile(int k, bool cmpl){
	for(int i=0; i<profWithFlanksLength; i++){
		Fourier *fx= (k==1) ? &kern->fy: &kern->fx;
		double  ReX=fx->re[i],
				ImX=fx->im[i],
				ReY=kern->ft.re[i],
				ImY=kern->ft.im[i];
		ImX=cmpl ? -ImX : ImX; ImY=cmpl ? -ImY : ImY;
		LCorrelation.datRe[i]=(ReX*ReY+ImX*ImY);
		LCorrelation.datIm[i]=(ReX*ImY-ImX*ReY);
	}
	LCorrelation.calc(0);				//reverse transformation
}
//==================== Make correlation track
double LocalCorrTrack(int pos, bool cmpl1, bool cmpl2){
	if(outWIG==NONE) return 0;
	getMem0(smoothProf2,profWithFlanksLength+10, "storeCorrTrack");

	calcSmoothProfile(1, cmpl2);	// calculate smooth ptrofile c=\int f*\rho for the second profile (y)
	memcpy(smoothProf2,LCorrelation.re,profWithFlanksLength*sizeof(double));
	calcSmoothProfile(0, cmpl1);	// calculate smooth profile for the first profile (x)
	smoothProf1=LCorrelation.re;

	double av=0;
	double sd=bTrack1.sd0*bTrack2.sd0;
	for(int i=LFlankProfSize; i<profWithFlanksLength-RFlankProfSize; i++){
		double x=smoothProf1[i]	/profWithFlanksLength;					//the smoothed profile for x
		double y=smoothProf2[i]	/profWithFlanksLength;					//the smoothed profile for y
		double lc=0;														//local correlation
		if((outWIG & WIG_CENTER) == WIG_CENTER) {x-=bTrack1.av0; y-=bTrack2.av0; }
		if((outWIG & WIG_MULT) == WIG_MULT){
			lc=x*y/sd;
		}
		if((outWIG & WIG_SUM ) == WIG_SUM ){
			double xx=kern->fx.datRe[i], yy=kern->fy.datRe[i];
			if((outWIG & WIG_CENTER) == WIG_CENTER) {xx-=bTrack1.av0; yy-=bTrack2.av0; }
			lc=0.5*(xx*y+yy*x)/sd;
		}
		LCorrelation.im[i]=lc; av+=lc;	// We use wCorrelation.re as a tmp buffer
		dHist.add(lc,pos < 0 ? 1:0);
	}
	av/=profileLength;
	if(pos>=0) wigCorr.addArray(LCorrelation.im+LFlankProfSize,pos);

	return av;
}

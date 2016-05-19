/*
 * outwig.cpp
 * Module does output of the local correlations into wig file.
 *
 *  Created on: Dec 11, 2013
 *      Author: mironov
 */
#include "track_util.h"

struct ScaledAray{
	/**
	 * The class encodes double values into short int
	 */
	float *array;		// encoded array
	char  *counts;		// counts of observations in given points
	int curPos;			// the lsat position of filled part of the array
	double min,max;		// min and max values to provide encoding and decoding
	ScaledAray();
	void init();		// initiation : allocate array. Size = (genome size) / stepSize
	void addArray(double *f, int pos);		// add array of doubles to the scaled array
//	void rescale(double amin, double amax);	//
//	short scale(double x);					// convert double to short int using current min and max
//	short scale(double x, double amin, double amax);  // convert double to short using given min and max
//	double unscale(short k);				// convert given short to double using current min and max
	void write();							// write output files
};

ScaledAray wigCorr;
double *corr2;


//==========================================================================
ScaledAray::ScaledAray(){
	min=1.e+128; max=-1.e+128; curPos=0; array=0; counts=0;
}

void ScaledAray::init(){
	getMem0(array,profileLength+10, "ScaledAray::init #1");
	getMem0(counts,profileLength+10, "ScaledAray::init #2");
	zeroMem(counts,profileLength+10);
}


void ScaledAray::addArray(double *f, int pos){

	for(int i=0, j=pos; i<wProfSize; i++,j++){
		double x=f[i];
		if(j<curPos){
			x=(array[j]*counts[j]+x)/(counts[j]+1);
		}
		array[j]=x;
		counts[j]++;
	}
	curPos=pos+wProfSize;
}

void ScaledAray::write(){
	char bf[1024];
	sprintf(bf,"%s.wig",outFile);
	FILE *outWigFile=xopen(bf,"wt");
	verb("\nwrite WIG correlations");
	ScoredRange pos, pos0;

	min=1.e+128; max=-1.e+128;
	for(int i=0; i<profileLength; i++){
		if(min > array[i]) min = array[i];
		if(max < array[i]) max = array[i];
	}

double lmax=-1000;
	fprintf(outWigFile,"track type=wiggle_0 ");
	char *s=strrchr(outFile,'/'); if(s==0) s=outFile; else s++;
	fprintf(outWigFile,"description=\"correlation: %s\" \n",s);
	fprintf(outWigFile,"#Param ");
	if((outWIG & WIG_BASE)  == WIG_BASE  ) fprintf(outWigFile,"BASE "  );
	if((outWIG & WIG_CENTER)== WIG_CENTER) fprintf(outWigFile,"CENTER ");
	if((outWIG & WIG_SUM)   == WIG_SUM   ) fprintf(outWigFile,"SUM "   );
	if((outWIG & WIG_MULT)  == WIG_MULT  ) fprintf(outWigFile,"MULT "  );
	fprintf(outWigFile,"\n");
	fprintf(outWigFile,"#Scale data: min=%.4f;  max=%.4f;   LOG scale: x=1000*log(a*array[i]+b); a=(2.718281828-1)/(max-min); b=1-min*a\n", min,max);
	fprintf(outWigFile,"#Source statstics: av1=%.4f;  av2=%.4f; sd1=%.4f; sd2=%.4f\n", bTrack1.av0,bTrack2.av0, bTrack1.sd0,bTrack2.sd0);
int kmin=1000, kmax=0;
	double a=(2.718281828-1)/(max-min), b=1-min*a;
	for(int i=0; i<profileLength; i++){
		double z=a*array[i]+b;
		if(lmax < z) lmax=z;
		int k=(int)(log(z)*1000);
		if(kmin > k) kmin=k;
		if(kmax < k) kmax=k;
		double v1=bTrack1.getValue(i,0), v2=bTrack2.getValue(i,0);
		if(k>=outThreshold && v1*v2!=0){
			filePos2Pos(i,&pos,stepSize);
			if(pos.chrom != pos0.chrom  || pos.beg-pos0.beg > stepSize){
				fprintf(outWigFile,"fixedStep chrom=%s start=%li step=%i span=%i\n",
						pos.chrom,pos.beg,stepSize,stepSize);
			}
			pos0.chrom=pos.chrom; pos0.beg=pos.beg;
			fprintf(outWigFile,"%i\n",k);
		}
	}
	fclose(outWigFile);
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
		wCorrelation.datRe[i]=(ReX*ReY+ImX*ImY);
		wCorrelation.datIm[i]=(ReX*ImY-ImX*ReY);
	}
	wCorrelation.calc(0);
}
//==================== Make correlation track
double storeCorrTrack(int pos, bool cmpl1, bool cmpl2){
	if(outWIG==NONE) return 0;
	getMem0(corr2,profWithFlanksLength+10, "storeCorrTrack");

	kern->fx.re[0]=kern->fx.re0;
	kern->fy.re[0]=kern->fy.re0;
	kern->ft.re[0]=kern->ft.re0;

	calcSmoothProfile(1, cmpl2);
	memcpy(corr2,wCorrelation.re,profWithFlanksLength*sizeof(double));
	calcSmoothProfile(0, cmpl1);

	double av=0;
	double sd=bTrack1.sd0*bTrack2.sd0;
	for(int i=LFlankProfSize; i<profWithFlanksLength-RFlankProfSize; i++){
		double x=wCorrelation.re[i]/profWithFlanksLength;
		double y=corr2[i]          /profWithFlanksLength;
		if((outWIG & WIG_CENTER) == WIG_CENTER) {x-=bTrack1.av0; y-=bTrack2.av0;}
		if((outWIG & WIG_MULT) == WIG_MULT) wCorrelation.re[i]=x*y/sd;
		if((outWIG & WIG_SUM ) == WIG_SUM ){
			double xx=kern->fx.datRe[i], yy=kern->fy.datRe[i];
			if((outWIG & WIG_CENTER) == WIG_CENTER) {xx-=bTrack1.av0; yy-=bTrack2.av0;}
			wCorrelation.re[i]=0.5*(xx*x+yy*y)/sd;
		av+=wCorrelation.re[i];
		}
	}
	av/=profileLength;
	wigCorr.addArray(wCorrelation.re,pos);
	return av;
}

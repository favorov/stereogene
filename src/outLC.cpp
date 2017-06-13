/*
 * outwig.cpp
 * Module does output of the local correlations into wig file.
 *
 *  Created on: Dec 11, 2013
 *      Author: mironov
 */
#include "track_util.h"

double *smoothProf2, *smoothProf1;
FloatArray *lcProfile=0;
float L_lcTreshold,R_lcTreshold;

struct ProfileHist:DinHistogram{
	double *lFDR;		// left lFDR
	double *rFDR;		// right lFDR
	double *rCDF_obs; 	// 1-CDF
	double *rCDF_exp;
	double *lCDF_obs;	// CDF
	double *lCDF_exp;

	ProfileHist(int l);
	~ProfileHist();
	void print(FILE*f);
	void fin();
};

ProfileHist dHist(100000);
ProfileHist NormWHist(500);

double LClogScale=1000;

float normLC(float lc){	//normalize Local correlation with LCScale parameter
	if(LCScale==LOG_SCALE)
		return lc >=0 ? log(1+lc) : -log(1-lc);
	return lc;
}

void renormDistrib(){	// Produce a NormWHist distribution using dHist
	for(int i=0; i<dHist.l; i++){
		int counts[2]={dHist.cnts[0][i],dHist.cnts[1][i]};
		double v=dHist.getValue(i);
		for(int tt=0; tt<2; tt++){
			if(counts[tt]){
				NormWHist.add(v,counts[tt],tt);
			}
		}
	}
	NormWHist.fin();
	//================================== Calculate thresholds ==========================
	L_lcTreshold=NormWHist.getValue(0); R_lcTreshold=NormWHist.getValue(NormWHist.l-1);
	for(int i=0; i<NormWHist.l; i++){
		if(NormWHist.rFDR[i] >= lcFDR)  R_lcTreshold=NormWHist.getValue(i);
		else break;
	}
	for(int i=NormWHist.l-1; i>=0; i--){
		if(NormWHist.lFDR[i] >= lcFDR)  L_lcTreshold=NormWHist.getValue(i);
		else break;
	}
}


ProfileHist::ProfileHist(int l):DinHistogram(l){//========= profile histogrgamm initiation
	getMem(lFDR,l,"WigHist:getFDR");
	getMem(lCDF_obs,l,"WigHist:getFDR");
	getMem(lCDF_exp,l,"WigHist:getFDR");
	getMem(rFDR,l,"WigHist:getFDR");
	getMem(rCDF_obs,l,"WigHist:getFDR");
	getMem(rCDF_exp,l,"WigHist:getFDR");
}
ProfileHist::~ProfileHist(){
	xfree(lFDR,"WigHist:getFDR");
	xfree(lCDF_obs,"WigHist:getFDR");
	xfree(lCDF_exp,"WigHist:getFDR");
	xfree(rFDR,"WigHist:getFDR");
	xfree(rCDF_obs,"WigHist:getFDR");
	xfree(rCDF_exp,"WigHist:getFDR");

}

void ProfileHist::fin(){
	DinHistogram::fin();					// normalize results
	double obs=0, exp=0;
	for(int i=l-1; i>=0; i--){
		rCDF_obs[i]=(obs+=hist[0][i]*bin);  // observed
		rCDF_exp[i]=(exp+=hist[1][i]*bin);	// expected
		double psi=1/n[0]*bin;				// psudocounts
		double x=(exp+psi)/(obs+psi);
		rFDR[i]=x>1 ? 1:x;					// right FDR
	}
	obs=exp=0;
	for(int i=0; i<l; i++){
		lCDF_obs[i]=(obs+=hist[0][i]*bin);  // observed
		lCDF_exp[i]=(exp+=hist[1][i]*bin);	// expected
		double psi=1/n[0]*bin;				// psudocounts
		double x=(exp+psi)/(obs+psi);
		lFDR[i]=x>1 ? 1:x;					// left FDR
	}
}


void ProfileHist::print(FILE* f){						// print the histogram
	fprintf(f,"#  min=%.3f max=%.3f \n",min,max);
	fprintf(f,"#  hMin=%.2ef  hMax=%.2ef bin=%.3f\n",hMin,hMax,bin);
	fprintf(f,"#  Observed: e0=%.3f sd0=%.3f n0=%i\n",e[0],sd[0],n[0]);
	fprintf(f,"#  Expected: e1=%.3f sd1=%.3f n1=%i\n",e[1],sd[1],n[1]);
	fprintf(f,"value\tobs\tnObs\texp\tnExp\tr_CDF_obs\tr_CDF_exp\tR_FDR");
	fprintf(f,"\tl_CDF_obs\tl_CDF_exp\tL_FDR\n");
	int i0=0, i1=l;
	for(int i=0; i<l; i++){
		if(i1==l && (hist[0][i]!=0 || hist[1][i]!=0)) i0=i;
		if(hist[0][i]!=0 || hist[1][i]!=0) i1=i;
	}

	for(int i=i0; i<=i1; i++){
		double h0=hist[0][i];
		double h1=hist[1][i];
		double fdrR=rFDR[i]*100, fdrL=lFDR[i]*100;
		fprintf(f,"%.3f\t%.2e\t%6i\t%.2e\t%6i\t%.2e\t%.2e\t%.2f", getValue(i),h0,cnts[0][i],h1,cnts[1][i],
				rCDF_obs[i], rCDF_exp[i],fdrR);
		fprintf(f,"\t%.2e\t%.2e\t%.2f\n", lCDF_obs[i], lCDF_exp[i],fdrL);
	}
}

//==========================================================================
//=== put the data to the Local correlation profile.
//=== The method takes into account possible window overlapping
void addLCProf(double *f, int pos){
	int d=wProfSize-wProfStep;
	for(int i=0, j=pos; i<wProfSize && j < profileLength; i++,j++){
		float x=f[i];
		if(pos > 0 && i < d) x/=2;
		if(pos+wSize < profileLength && (i > wStep)) x/=2;
		lcProfile->add(j,x);
	}
}


//===================== Write the local correlation into the bedGraph file =====
void writeLC(){
	char bf[1024];
	verb("\nwrite Local Correlations\n");
	renormDistrib();
	writeBedGr(outFile, lcProfile, L_lcTreshold,R_lcTreshold);

	//========================================== Write histograms ===============
	if(writeDistr) {
		sprintf(bf,"%s.LChist",outFile);
		FILE *lcHistF=fopen(bf,"wt");
		NormWHist.print(lcHistF);
		fclose(lcHistF);
	}
	dHist.clear();
	NormWHist.clear();
}
//==================================================================
void initOutLC(){
	if(outLC==NONE) return;
	if(lcProfile==0) lcProfile=new FloatArray();
	lcProfile->init(NA);
}

void finOutLC(){
	if(outLC==NONE) return;
	//====================== write correlation
	writeLC();
}

void freeLC(){
	if(lcProfile) delete lcProfile;
	lcProfile=0;
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
		LCorrelation.datRe[i]= (ReX*ReY+ImX*ImY);
		LCorrelation.datIm[i]= (ReX*ImY-ImX*ReY);
	}
	LCorrelation.datRe[0]=0;
	LCorrelation.datIm[0]=0;

	LCorrelation.calc(0);				//reverse transformation
}
//==================== Make correlation track
double *lcTmp=0;

double LocalCorrTrack(int pos1, int pos2, bool cmpl1, bool cmpl2, bool rnd){

	if(outLC==NONE) return 0;
	getMem(smoothProf1,profWithFlanksLength+10, "storeCorrTrack");
	getMem(lcTmp,profWithFlanksLength+10, "storeCorrTrack");
	calcSmoothProfile(0, cmpl1);	// calculate smooth ptrofile c=\int f*\rho for the second profile (y)
	memcpy(smoothProf1,LCorrelation.re,profWithFlanksLength*sizeof(double));
	calcSmoothProfile(1, cmpl2);	// calculate smooth profile for the first profile (x)
	smoothProf2=LCorrelation.re;

	double av=0;
	double sd=track1->sd0*track2->sd0;
	for(int i=LFlankProfSize; i<profWithFlanksLength-RFlankProfSize; i++){
		double x=smoothProf1[i]	/profWithFlanksLength;					//the smoothed profile for x
		double y=smoothProf2[i]	/profWithFlanksLength;					//the smoothed profile for y
		double lc=0;													//local correlation
//		int ax=bTrack1.getBVal(pos1+i-LFlankProfSize,cmpl1);
//		int ay=bTrack2.getBVal(pos2+i-LFlankProfSize,cmpl2);
//		if(ax==NA || ay==NA) {lcTmp[i]=NA; continue;}
//		double x0=(ax==NA)?0: bTrack1.getVal(ax), y0=(ay==NA)?0: bTrack2.getVal(ay);
//		lc=0.5*(x0*y+y0*x);
		lc=x*y/sd;
//if(!rnd) deb("%i  bb=%i,%i vxy=%f,%f xy=%f,%f lc=%f",i,ax,ay,bTrack1.getVal(ax), bTrack2.getVal(ay),x,y,lc);
		lc=normLC(lc);
		lcTmp[i]=lc; av+=lc;	    // We use wCorrelation.im as a tmp buffer
		dHist.add(lc,rnd ? 1:0);
	}
	av/=profileLength;
	if(!rnd) addLCProf(lcTmp+LFlankProfSize,pos1);
	return av;
}

//====== write bed graph file. lTreshold,rTreshold -- thresholds (NA allowed).
//       The value will be written if it is OUTSIDE the threshold interval
void writeBedGr(FILE* f, FloatArray *array){
	writeBedGr(f,array,NA,NA);
}

void writeBedGr(const char *fname, FloatArray *array, float lTreshold, float rTreshold){
	char b[1024];
	makeFileName(b,resPath,fname);
	strcat(strcat(b,"."),BGR_EXT);
	FILE *f=fopen(b,"w");
	char bname[1024];
	strcpy(bname,outFile);
	char *s=strrchr(bname,'/'); s= s==0 ? bname : s+1;
	fprintf(f,"track type=bedGraph name=\"%s\" description=\"Local correlation. FDR=%.1f%%\"\n", s, lcFDR*100);
	writeBedGr(f,array, lTreshold,  rTreshold);
	fclose(f);
}

void writeBedGr(FILE* f, FloatArray *array, float lTreshold, float rTreshold){
	ScoredRange pos, pos0;
	if(lTreshold == NA) lTreshold= 1.e+8;
	if(rTreshold == NA) rTreshold=-1.e+8;
	for(int i=0; i<profileLength; i++){
		if(i%1000000 ==0) verb("%5.1f%%\r",1.*i/profileLength*100);
		float x=array->get(i);
		if(x <= rTreshold && x >= lTreshold) x=NA;
		filePos2Pos(i,&pos,binSize); pos.score=x;
		if(pos.chrom == pos0.chrom && pos.score == pos0.score){
			pos0.end=pos.end; continue;
		}
		if(pos0.chrom !=0)
			pos0.printBGraph(f);
		pos0.beg=pos.beg; pos0.end=pos.end; pos0.chrom=pos.chrom; pos0.score=pos.score;

	}
	if(pos0.chrom !=0) pos0.printBGraph(f);
}



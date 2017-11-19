/*
 * main_proj.cpp
 *
 *  Created on: 03 ����. 2017 �.
 *      Author: andrey
 */
#include "track_util.h"

const char * progName="Projection";
const int progType=SG;


void printProgDescr(){
	printf("\n");
	printf("The Smoother program creates smoothed track using a kernel\n");
	printf("Usage:\n");
	printf("$ ./Smoother [-parameters] track1 track2 ...\n");
	printf("\n");
}
void printMiniHelp(){
	printf("\n");
	printf("The Smoother program creates smoothed track using a kernel\n");
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./Smoother [-parameters] track1 track2 ...\n");
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}
//=====================================================================
double *smoothProf=0;
double *smTmp=0;


void smooth(const char *fname){
	outLC=1;
	initOutLC();
	int l=profileLength;
	bTrack *tr=new bTrack(fname);

	for(int i=0,k=0; i<l; i+=wProfStep,k++){
		double d;
		d=100.*k/(l/wProfStep);
		if(k%10000 ==0) verb("\nSmooother: %4.1f%% (%6i/%i) ",d,k,l/wProfStep);
		else if(k%1000 ==0) verb(".");
		double *pr1=tr->getProfile(i,0);		// decode the first profile. Decoder uses hasCompl and complFg flags and combines profiles
		if(smoothProf==0) getMem(smoothProf,profWithFlanksLength+10, "storeCorrTrack");
		if(smTmp==0)      getMem(smTmp,profWithFlanksLength+10, "storeCorrTrack");
		kern->fftx(pr1,0);
		calcSmoothProfile(&(kern->fx),0, 0);	// calculate smooth ptrofile c=\int f*\rho
		addLCProf(LCorrelation.re+LFlankProfSize,i);
	}
	char pfil[4096],wfil[4096];
	makeFileName(pfil,trackPath,fname);

	char *s=strrchr(pfil,'/'); if(s==0) s=wfil;
	s=strrchr(s,'.'); if(s) *s=0;
	sprintf(wfil,"%s_sm.bgr",pfil);

	//================ normalize
	double tt=0;
	for(int i=0; i<l; i++) tt+=lcProfile->get(i);
	for(int i=0; i<l; i++) {
		double x=lcProfile->get(i)*tr->total/tt;
//		if(x<0.5) x=0;
		int xx=x*50;
		x=xx/50.;
		lcProfile->set(i,x);
	}


	FILE *f=gopen(wfil,"w");
	fprintf(f,"track type=bedGraph name=\"%s_Smooth\" description=\"Smoothed track. Width=%.0f\"\n",fname, kernelSigma);
	verb("\n Write Smooth profile...\n");
	writeBedGr(f,lcProfile, NA,  NA);
	verb("   Done\n");
}

//=====================================================================
void Smoother(){
	PrepareParams();
	LCorrelation.init(profWithFlanksLength);
	lcFlag=BASE;
	getMem(LCorrelation.datRe,profWithFlanksLength, "Correlator");


	for(int i=0; i<nfiles; i++){
		char *fname=files[i].fname;
		if(fname==0 || strlen(trim(fname))==0) continue;
		smooth(fname);
	}
}

//=====================================================================
int main(int argc, char **argv) {
	initSG(argc, argv);
	if(debugFg) {clearDeb(); debugFg=DEBUG_LOG|DEBUG_PRINT;}

	Preparator();

	Smoother();
	fflush(stdout);
	fclose(stdout);
	return 0;
}





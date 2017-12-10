/*
 * data.cpp

 *
 *  Created on: 05 дек. 2017 г.
 *      Author: andrey
 */

#include "track_util.h"

const char* version="2.13";


Chromosome *chrom_list;       // list of chromosomes
Chromosome *curChrom=chrom_list;
int  binSize=100;   // frame size fo profile
bool  NAFlag=0;

long long GenomeLength=0;      // TOTAL LENGTH OF THE GENOME
int n_chrom;
char trackName[4096];       // current track name
char *chromFile=0;
char *confFile=0;		// confounder file name
char *cfgFile=0;		// config file name
char *profPath=0;
char *trackPath=0;

char *trackFil=0;		// Track file
char *profile1=0;		// first profile file file name
char *profile2=0;		// second profile file file name
char *resPath=0;
char *statFileName=(char*)"./statistics";
char *paramsFileName=(char*)"./params";
char *inputProfiles=0;
char *outTrackFile=0; // Filename for write out track
char *idSuff=(char*)"";

bool  syntax=1;				// Strong syntax control

bool  writeDistr=1;
bool  writeDistCorr=1;		    // write BroadPeak
int   crossWidth=10000;
bool  outSpectr=0;
bool  outChrom=0;
int   outRes=XML|TAB;
int   inpThreshold=0;		// Testing of binarized input data, % of max
bool  writePDF=true;
int   complFg=IGNORE_STRAND;
int   lcFlag=CENTER;
int   profileLength;			// size of the profile array

char *pcorProfile=0;    // partial correlation profile file name
int  binBufSize=30000000;


int 	kernelType=KERN_NORM;
char* 	customKern=0;
double 	noiseLevel=0.2;
int 	wSize=100000;        // size of widow (nucleotides)
int 	wStep=0;             // window step   (nucleotides)
int 	flankSize=0;
double 	kernelSigma=1000.;    	// kernel width (nucleotides)
double 	kernelShift=0;      	    // Kernel mean (for Gauss) or Kernel start for exponent
int 	intervFg0;
double 	scaleFactor=0.2;
bool 	outLC=0;
int 	LCScale=LOG_SCALE;
//int 	LCScale=LIN_SCALE;
double LlcFDR=0;		// treshold on FDR when write Local Correlation track
double RlcFDR=0.5;		// treshold on FDR when write Local Correlation track

bool 	outPrjBGr=true;

int 	wProfStep=0;          	// window step   (profile scale)
int 	wProfSize=0;          	// size of widow (profile scale)
int 	LFlankProfSize=0;         // size of flank (profile scale)
int 	RFlankProfSize=0;         // size of flank (profile scale)
int 	profWithFlanksLength=0; 	// size of profWindow array (including random flanks)
double 	kernelProfSigma=1000;     // kernel width ((profile scale)
double 	kernelProfShift=0;
double 	kernelNS=0;			// Correction for non-specifisity
Track 	*track1=0, *track2=0, *projTrack=0;
Kernel 	*kern=0;
double 	maxNA0=95;
double 	maxZero0=95;
double 	maxNA=100;
double 	maxZero=100;
int 	nShuffle=10000;
char *trackName1=strdup("");
char *trackName2=strdup("");
double mannW_Z=0;
double mannW_p=1;

Model 	*model;

int 	threshold=0;

FILE 	*logFile=0;
bool 	doAutoCorr=0;

int 	corrScale=10;
double 	prod11=0,prod12=0,prod22=0, eprod1,eprod2;
int 	nprod=0;
XYCorrelation XYfgCorrelation;		    // array for correlation picture
XYCorrelation XYbgcorrelation;			// array for correlation picture
Fourier LCorrelation;

double 	totCorr=0, BgTotal=0;
bool 	RScriptFg=0;
int 	bpType=BP_SIGNAL;
int 	cage=0;
bool 	clearProfile=false;
int 	scoreType=AV_SCORE;
FileListEntry files[256];
int   	nfiles;
bool LCExists=false;

double 	BgAvCorr=0;
double 	FgAvCorr=0;
int  	pgLevel=2;
float total=0;						// total count over the track

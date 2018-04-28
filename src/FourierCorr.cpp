/*
 * FourierCorr.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: mironov
 */
#include "track_util.h"

double *BkgSet=0, *FgSet=0;			// background and foreground sets of the correlations
int nBkg, nFg;					// size of background and foreground sets of data
int maxPairs=0, fstep=0;
double *xDat,*yDat,*xyCorr;  	//Woriking arrays

PairEntry *pairs;				// array for pair's correlation (foreground)
int nPairs;						// number of foreground observations

void clear(){
	writeLog("clear Correlator\n");
	if(FgSet) xfree(FgSet,"FgSet"); FgSet=0;
	if(pairs) xfree(pairs,"pairs"); pairs=0;
	nBkg=0; nFg=0; nPairs=0;
}

//============================================== subtract profile from profile
double scalar(double *px, double *py, int l){	//partial correlation variant 1
	double res = 0;
	for (int i = 0; i < l; i++){
		res += px[i]*py[i];
	}
	return res;
}

void minusProf(double *profx, double *profz, double pcorCoef){	//partial correlation variant 1

	for (int i = LFlankProfSize; i < LFlankProfSize + profileLength; i++){
		profx[i] -= profz[i]*pcorCoef;
	}
}


//============================================== Calculate correlation for given pair of windows
//============================================== Add statistics to chromosome
void addChromStat(int pos, double corr, double lCorr, double av1, double av2){
	ScoredRange gr;

	filePos2Pos(pos,&gr,0);

	Chromosome *chr=gr.chr;
	chr->corr+=corr;
	chr->lCorr+=lCorr;
	chr->av1+=av1;
	chr->av2+=av2;
	chr->count++;
}
//=======  input: to position and complement flags. rnd: if the data comes from shuffling

double calcCorelations(int pos1, int pos2, bool cmpl1, bool cmpl2, bool rnd){
	int na1=track1->countNA(pos1,cmpl1);			// count Na's  in the first profile
	int na2=track2->countNA(pos2,cmpl2);			// count Na's  in the second profile
	int nz1=track1->countZero(pos1,cmpl1);			// count zeros in the first profile
	int nz2=track2->countZero(pos2,cmpl2);			// count zeros in the second profile
	if(na1 > maxNA) {return -101;}					// too many NA in the first profile
	if(na2 > maxNA) {return -102;}					// too many NA in the second profile
	if(nz1 > maxZero) return -201; 					// too many zeros in the profiles
	if(nz2 > maxZero) return -202; 					// too many zeros in the profiles

	double *pr1=track1->getProfile(pos1,cmpl1);		// decode the first profile. Decoder uses hasCompl and complFg flags and combines profiles
	double *pr2=track2->getProfile(pos2,cmpl2);		// decode the second profile
	kern->fftx(pr1,track1->deriv);					// do fft for the profiles
	kern->ffty(pr2,track2->deriv);
	double corr=kern->dist(cmpl1);					// Kernel strand is selected by the first profile
	double lCorr=0, av1, av2;
	if(corr > -10){									// Error in the correlation => skip the pair of the windows
		if(outLC){				                 // Make local correlation track
			lCorr=LocalCorrTrack(pos1, pos2, cmpl1, cmpl2,rnd);
		}
		if(!rnd) {									// foreground distribution
			//======= Calc the correlation
			if(writeDistCorr) XYfgCorrelation.calcXYCorr(pos1,cmpl1, cmpl2,corr);
			if(outSpectr    ) XYfgCorrelation.makeSpectrum();

			av1=track1->addStatistics();
			av2=track2->addStatistics();
			addChromStat(pos1,corr,lCorr, av1,av2);	// Add data to chromosome statistics
			if(doAutoCorr){calcAutoCorr();}
		}
		else{
			if(writeDistCorr)
				{
				XYbgcorrelation.calcXYCorr(-1,cmpl1, cmpl2,corr); //Do background correlation function
				}
		}
	}
	return corr;
}
//======================================================= Calculate the total correlation
double calcCC(){
	double c11=prod11-eprod1*eprod1/nprod*profWithFlanksLength;
	double c22=prod22-eprod2*eprod2/nprod*profWithFlanksLength;
	double c12=prod12-eprod1*eprod2/nprod*profWithFlanksLength;
	double cc=c12/sqrt(c11*c22);
	return cc;
}

void cleanCummulative(){
	prod11=0; prod12=0; prod22=0; eprod1=0; eprod2=0; nprod=0;
}
//======================================================== Calculate background distributions
struct IntPair{
	int p1,p2;
	char cmpl1, cmpl2;
	IntPair(){
		cmpl1=track1->hasCompl && (rand() > RAND_MAX/2);
		if(cmpl1){if(track1->ivsC->nIv==0) cmpl1=0;}
		else     {if(track1->ivs ->nIv==0) cmpl1=1;}

		cmpl2=track2->hasCompl && (rand() > RAND_MAX/2);
		if(cmpl2){if(track2->ivsC->nIv==0) cmpl2=0;}
		else     {if(track2->ivs ->nIv==0) cmpl2=1;}
		p1=track1->getRnd(cmpl1);
		p2=track2->getRnd(cmpl2);
	}
};

IntPair *posPairs;

int posPairCmp(const void *xp1, const void *xp2){
	IntPair *pair1=(IntPair *) xp1;
	IntPair *pair2=(IntPair *) xp2;
	if(pair1->cmpl1 != pair2->cmpl1) return (int)(pair1->cmpl1) - (int)(pair2->cmpl1);
	if(pair1->cmpl2 != pair2->cmpl2) return (int)(pair1->cmpl2) - (int)(pair2->cmpl2);

	int bl1=pair1->p1/binBufSize;
	int bl2=pair2->p1/binBufSize;
	if(bl1!=bl2) return bl1-bl2;
	return pair1->p2 -pair2->p2;
}

int distrBkg(int n);

void distrBkg(){
	srand(33);									// random seed
	verb("\nBakcground...");
	cleanCummulative();
	avBg=0;
	getMem0(posPairs,nShuffle,"init randomPairs");
	getMem0(BkgSet,nShuffle, "bkg Distr"); nBkg=0; 	// allocate array for background observations
	int n=nShuffle;
	int tst=0;
	do{
		distrBkg(n); n=nShuffle-nBkg;
		if(tst++ > 100) {
			writeLog("Too few non-zero windows. Only %i shuffles done", nBkg);
			break;
		}
	}while(n);

	double cc=calcCC();
	avBg/=nBkg;
	BgTotal=cc;
	xverb("\nbg_cc=%f \nbg_average=%f\n",cc,avBg);
	errStatus=0;
}


int distrBkg(int nSh){
	errStatus="bkg Distrib.";
	//================================================== simulations
	int tst=0;											// count runs with zero/NA frames
	for(int i=0; i<nSh; i++)	posPairs[i]=IntPair();
	qsort(posPairs,nSh,sizeof(IntPair),posPairCmp);
	for(int i=0; i<nSh; i++){
		int  p1   =posPairs[i].p1, p2=posPairs[i].p2;
		bool cmpl1=posPairs[i].cmpl1, cmpl2=posPairs[i].cmpl2;

		double d=calcCorelations(p1,p2, cmpl1, cmpl2, true);		// calculate correlation
		if(d<=-10) {							// invalid windows (too many NA's of Zeros)
			if(tst++ > 10000){					// too many attempt to get a background correlations
				writeLog("too many empty/zero windows\n"); return tst;
			}
			continue;
		}
		avBg+=d;
		if(i%1000 ==0) verb("\nShuffling: %i/%li",i,nShuffle);
		else if(i%100 ==0) verb(".");
		tst=0;
		BkgSet[nBkg++]=d;						// store in distribution
	}
	return tst;
}


//============================================ Store foreground distribution
inline void storePair(int i, double d){
	FgSet[nFg++]=d;
	if(writeDistr==DISTR_DETAIL){
		PairEntry *pe=pairs+(nPairs++);						//== store pair of positions
		pe->profPos=i; pe->d=(float)d;						//== define the pair
	}
}

//============================================= Calculate coherent correlations
void distrCorr(){
	verb("\nForeground...");
	int l=profileLength;
	errStatus="distrCorr";
	maxPairs=l/wProfStep; fstep=wProfSize/wProfStep;
	if(fstep==0) fstep=1;

	int nTrkPair=1; if(track1->hasCompl) nTrkPair*=2; if(track2->hasCompl) nTrkPair*=2;
	maxPairs*=nTrkPair;

	int siz=(maxPairs+100);
	getMem0(FgSet, siz, "Corr #1");	zeroMem(FgSet, siz);		//== array for foreground distribution
	if(writeDistr==DISTR_DETAIL) {getMem0(pairs, siz, "Corr #2");	zeroMem(pairs, siz);}		//== array for pairs
	cleanCummulative();

	//=================== calculate correlations
	int n_corr=0; avFg=0;

	for(int i=0,k=0; i<l; i+=wProfStep,k++){
		double d;
		d=100.*k/(l/wProfStep);
		if(k%10000 ==0) verb("\ncoherent: %4.1f%% (%6i/%i) ",d,k,l/wProfStep);
		else if(k%1000 ==0) verb(".");

		if((complFg==IGNORE_STRAND)||(!track1->hasCompl && !track2->hasCompl)){ // no direction defined
			if((d=calcCorelations(i,i, false,false,false)) >=-10){
				storePair(i,d); n_corr++; avFg+=d;
			}
		}
		else if((complFg&COLLINEAR)!=0){								// analyze collinear chains
			if((d=calcCorelations(i,i,true,true,false)) >=-10){		// => =>  valid pair
				storePair(i,d); n_corr++; avFg+=d;
			}
			if((d=calcCorelations(i,i,false,false,false)) >=-10){	// <= <=  valid pair
				storePair(i,d); n_corr++; avFg+=d;
			}
		}
		else if((complFg&COMPLEMENT)!=0){						// analyze complement chains
			if((d=calcCorelations(i,i,true,false,false)) >=-10){		// => =>  valid pair
				storePair(i,d); n_corr++; avFg+=d;
			}
			if((d=calcCorelations(i,i,false,true,false)) >=-10){	// <= <=  valid pair
				storePair(i,d); n_corr++; avFg+=d;
			}
		}
	}					// end for

	//=================================================== Define rank for q-value calculation
	if(n_corr==0){
		xverb("\nno non-zero windows pairs\n",totCorr, avFg);
	}
	else{
		track1->finStatistics(); track2->finStatistics();
		avFg/=n_corr;
		finOutLC();
		totCorr=calcCC();
		xverb("\nCorrelation=%f\naverage Corrrelation=%f\n",totCorr, avFg);
	}
	errStatus=0;
}

void calcAutoCorr(){
	double *autoX=	kern->fx.getAutoCorr();
	double *autoY=	kern->fy.getAutoCorr();
	for(int i=0; i<profWithFlanksLength; i++){
		track1->autoCorr[i]+=autoX[i];
		track2->autoCorr[i]+=autoY[i];
	}
}

//================================================================================================
//================================================================================================
//================================================================================================
//================================================================================================
char *resFileName(const char* n1,const char* n2){
	char b[2048];
	sprintf(b,"%s~%s",n1,n2);
	return strdup(b);
}
//================================================================== Make name for outfile
char * makeOutFilename(char * prof1, char*prof2){
	char p1Fname[4096], p2Fname[4096], b[4096];
	getFnameWithoutExt(b, prof1);
	if(strchr(b,'~')) sprintf(p1Fname,"(%s)",b);
	else strcpy(p1Fname,b);
	getFnameWithoutExt(b, prof2);
	if(strchr(b,'~')) sprintf(p2Fname,"(%s)",b);
	else strcpy(p2Fname,b);

	sprintf(b,"%s%s",resPath,p1Fname);
	return resFileName(b,p2Fname);
}

//==================== Check for duplication of the comparison
struct FilePair{
	FileListEntry *fil1, *fil2;
	FilePair(FileListEntry *f1, FileListEntry *f2){
		fil1=f1; fil2=f2;
	}
};
const int maxFilePairs=0x10000;
FilePair *fPairs[maxFilePairs];
int nFPairs=0;

int addPair(FileListEntry *f1, FileListEntry *f2){
	FileListEntry *ff1=f1,*ff2=f2;

	for(int i=0; i<nFPairs; i++){
		if(strcmp(ff1->fname, fPairs[i]->fil1->fname)==0 && strcmp(ff2->fname, fPairs[i]->fil2->fname)==0) {return 0;}
		if(strcmp(ff2->fname, fPairs[i]->fil1->fname)==0 && strcmp(ff1->fname, fPairs[i]->fil2->fname)==0) {return 0;}
	}
	if(nFPairs >= maxFilePairs) errorExit("too many comparisons in one run");
	fPairs[nFPairs++]=new FilePair(ff1,ff2);
	return nFPairs;
}

//========================================================================================
int fPaircmp(const void *a1, const void* a2){
	FilePair *p1=(FilePair *)a1;
	FilePair *p2=(FilePair *)a2;
	int x=strcmp(p1->fil2->fname, p2->fil2->fname);
	if(x) return x;
	return strcmp(p1->fil1->fname, p2->fil1->fname);
}

int Correlator(){
	Timer timer;
	id=0;	// id is undefined yet
	//================================================================== print parameters
	verb("========== Parameters ===========\n");
	if (pcorProfile != 0) verb("==         pcorProfile=<%s>\n", pcorProfile);
	verb("===        bin=%i\n",binSize);
	verb("==         wSize=%i\n",wSize);
	verb("==         kernelSigma=%.0f\n",kernelSigma);
	verb("==         nShuffle=%i\n",nShuffle);

	PrepareParams();
	//============ Read Map File

	if(pcorProfile) {
		projTrack=new bTrack();
		verb("read confounder...\n");
		projTrack->openTrack(pcorProfile);
	}
	int n_cmp=0;
	int nnf=nfiles; if(nnf>1) nnf--;
	LCorrelation.init(profWithFlanksLength);
	getMem0(LCorrelation.datRe,profWithFlanksLength, "Correlator");

	//============================================= Make Profile Pairs
	for(int i=0; i< nfiles; i++){
		for(int j=i+1; j< nfiles; j++){
			if(files[j].listId==files[i].listId) continue;
			if(strcmp (files[j].fname,files[i].fname)==0) continue;
			if(addPair(files+i,files+j)==0) continue;
		}
	}

	qsort(fPairs, nPairs, sizeof(FilePair), fPaircmp);
	//============================================== Do comparison
	FileListEntry *fil1=0, *fil2=0;

	for(int i=0; i<nFPairs; i++){
		id=(unsigned int)time(0);					// id is undefined yet

		if(fPairs[i]->fil1 != fil1){
			fil1=fPairs[i]->fil1;
			trackName1=fil1->fname;
			verb("read profile1 <%s>\n", trackName1);
			if(track1) {del(track1); track1=0;}
			track1=trackFactory(trackName1);
			if(pcorProfile) track1->ortProject();
			if(!track1->makeIntervals()){continue;}
		}
		if(fPairs[i]->fil2 != fil2){
			fil2=fPairs[i]->fil2;
			trackName2=fil2->fname;
			verb("read profile2 <%s>\n", trackName2);
			if(track2) {del(track2); track2=0;}
			track2=trackFactory(trackName2);
			if(pcorProfile) track2->ortProject();
			if(!track2->makeIntervals()) {continue;}
		}
		Timer thisTimer;
		outFile=makeOutFilename(trackName1, trackName2);
		makeId();
		writeLog("Correlations:   wSize=%i   kernelType=\"%s\"   kernelSigma=%.0f\n",
				wSize,getKernelType(),kernelSigma);
		writeLog("  in1=<%s> in2=<%s> out=<%s>\n", trackName1, trackName2, outFile);

		xverb("in1=\"%s\"\n", trackName1);
		xverb("in2=\"%s\"\n", trackName2);
		xverb("out=\"%s\"\n", outFile);
		//===================================================================== Calculate
		clearChromosomes();

		if(sparse){
			calcSparce();
		}
		else{
			XYfgCorrelation.initXY();
			XYbgcorrelation.initXY();
			initOutLC();
			writeLog("Background\n");
			distrBkg();						// Make background distribution

			writeLog("Foreground\n");
			distrCorr();					// Calculate correlations
			writeLog("Correlations -> Done\n");
		}
		if(nFg && nBkg){
			printCorrelations();			// write correlations
			printStat();					// write report
			if(RScriptFg) {
				printR();
				if(writePDF){
					printRreport();
					printRmd();
				}
			}
		}
		else{
			xverb("*** <%s>: No data for statistics: nFg=%i nBkg=%i ***\n", outFile,nFg, nBkg);
		}
		n_cmp++;
		clear();
		writeLog("<%s> => Done  time=%s\n",outFile,thisTimer.getTime());
	}
	freeLC();
	writeLog("====== DONE ======\n");
    verb("***   calculation time for %i comparisons = %s\n",n_cmp, timer.getTime());
	return 0;
}

/*
 * sparce.cpp
 *
 *  Created on: 27 февр. 2018 г.
 *      Author: andrey
 */
#include "track_util.h"

Track *master, *slave;  // master -- sparce track

void testDat(int l){
	for(int i=0; i<l; i++){
		deb("%i\t%.2e\t%.2e\t%.2e\t%.2e\t:\t%.2e\t%.2e\t%.2e\t%.2e\t:\t%.2e\t%.2e\t%.2e\t%.2e",i
				,kern->fx.datRe[i],kern->fx.datIm[i],kern->fx.re[i],kern->fx.im[i]
				,kern->fy.datRe[i],kern->fy.datIm[i],kern->fy.re[i],kern->fy.im[i]
				,kern->ft.datRe[i],kern->ft.datIm[i],kern->ft.re[i],kern->ft.im[i]
				);
	}
}


struct sparceDat{
	int pos1,pos2,l;
	double d11,d12,d22,e1,e2;
	sparceDat(int p1, int p2, int ll){pos1=p1; pos2=p2; l=ll; d11=d12=d22=e1=e2=0;}
	void calcDist();
};
//=====================================================================
int debbfg=0;
void sparceDat::calcDist(){
	int ll=l+LFlankProfSize+RFlankProfSize;
	master->getProfile(kern->fx.datRe, pos1, l, 0);
	slave->getProfile(kern->fy.datRe, pos2, l, 0);

	kern->fftx(kern->fx.datRe,0); kern->ffty(kern->fy.datRe,0);
if(debbfg)testDat(ll);
	kern->dist(0);
	d11=kern->dx11; d12=kern->dx12;	d22=kern->dx22; e1=kern->ex1; e2=kern->ex2;
if(debbfg) deb("%i\t%x\t%x\t%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e",pos1, 2*pos1, 2*pos2,l, d11, d12, d22, e1,e2);
}


//========================================================================
sparceDat** sparceCoherent;
sparceDat** sparceRand;

double calcSparceCorr(sparceDat** sp, int nn, int bootStrap){
	double e1=0,e2=0,d11=0,d12=0,d22=0,l=0;
	for(int i=0; i<nn; i++){
		if(bootStrap && drand() < 0.5) continue;
		sparceDat* sd=sp[i];
		l+=sd->l; e1+=sd->e1; e2+=sd->e2; d11+=sd->d11; d12+=sd->d12; d22+=sd->d22;
	}
	double cc=d12/sqrt(d11*d22);
//deb("==========l=%.0f e1=%.2e e2=%.2e d11=%.2e d12=%.2e d22=%.2e    cc=%f",l, e1, e2, d11, d12, d22,cc);
	return cc;
}

//=====================================================================================
int sparseCmp(const void* a1, const void *a2){
	sparceDat *d1=*((sparceDat**) a1);
	sparceDat *d2=*((sparceDat**) a1);
	return d1->pos2 - d2->pos2;
}

int calcSparce(){

	if(track1->ivs->totLength < track2->ivs->totLength){
		master=track1; slave=track2;
	}
	else{
		master=track2; slave=track1;
	}
	int nn=master->ivs->nIv;
	if(nn < 100) {errorExit("to few intrevals for sparce data analisys: %i\n", nn); return -1;}

	getMem(sparceCoherent, nn, "sparce1");
	getMem(sparceRand, nn, "sparce2");
	int lFlank=LFlankProfSize, rFlank=RFlankProfSize;	//save flank size

	//============================================== fill positions
	for(int i=0; i<nn; i++){
		int f=master->ivs->ivs[i]->f, t=master->ivs->ivs[i]->t;
		if(f < 0) f=0; if(t > profileLength) t=profileLength;
		int ll=t-f;
		sparceCoherent[i]=new sparceDat(f, f, ll);
		int rf=slave->getRnd(0);
		sparceRand[i]=new sparceDat(f,rf,ll);
	}
	qsort(sparceRand,nn,sizeof(sparceDat*),sparseCmp);
	//============================================== calc corr dat
//deb("p1\tp2\tl\td11\td12\td22\te1\te2");


	for(int i=0; i<nn; i++){
		if(i%1000==0) verb("%i / %i         \r",i,nn);
		sparceDat *sd=sparceCoherent[i];
		defFlanks(sd->l);
		int l=sd->l+LFlankProfSize+RFlankProfSize;
		getMem(kern->fx.datRe,l,"sparse kern#1"); getMem(kern->fy.datRe,l,"sparse kern#2");
		sparceCoherent[i]->calcDist();
		sparceRand    [i]->calcDist();
		xfree(kern->fy.datRe,0); xfree(kern->fx.datRe,0); kern->fx.datRe=kern->fy.datRe=0;
	}

	verb("\nBootstrap\n");
	//============================================ calc correlations
	nFg=nBkg=nShuffle;
	if(FgSet ==0) getMem(FgSet, nFg,"sparce #3");
	if(BkgSet==0) getMem(BkgSet, nBkg,"sparce #3");
	totCorr=calcSparceCorr(sparceCoherent,nn, 0);
	BgTotal=calcSparceCorr(sparceRand,nn, 0);
	for(int i=0; i<nShuffle; i++){		//======== Bootstrap
		if(i%100==0) verb("%i / %i         \r",i,nShuffle);
		double c1=calcSparceCorr(sparceCoherent,nn, 1);
		double c2=calcSparceCorr(sparceRand    ,nn, 1);
		FgSet [i]=c1;
		BkgSet[i]=c2;
	}
	xverb("\nCorrelation=%f\n",totCorr);

	//=================================================================================
	LFlankProfSize=lFlank; RFlankProfSize=rFlank;	//restore flank size
	for(int i=0; i<nn; i++) {delete sparceCoherent[i]; delete sparceRand[i];}
	free(sparceCoherent); sparceCoherent=0;
	free(sparceRand); sparceRand=0;
	return 0;
}




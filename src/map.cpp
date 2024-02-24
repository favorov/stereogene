/*
 * map.cpp
 *
 *  Created on: 09.03.2013
 *      Author: Mironov
 */
#include "track_util.h"






//===============================================================================
	MapRange::MapRange(){f=0; t=0; cumLength=0;}
	MapRange::MapRange(int ff, int tt){f=ff; t=tt; cumLength=0;}


//===============================================================================
IVSet::IVSet(){
	capacity=1000; nIv=0; totLength=0;
	errStatus="ivset init";
	getMem(ivs,capacity, "ivset init #1");
	errStatus=0;
}


IVSet::~IVSet(){
	for(int i=0; i<nIv; i++) delete ivs[i];
	xfree(ivs,"~ivs");
}




void IVSet::clear(){
	nIv=0; totLength=0;
	errStatus="ivset init";
	errStatus=0;
}
//===============================================================================
void IVSet::addIv(int f, int t){
	if(f < 0) {f = 0;}    if(t > profileLength) {t = profileLength;}
	if(t <= f) return;


	if(nIv > 0 && f < ivs[nIv-1]->t){//================ Overlaps: collect with previous
		if(ivs[nIv-1]->f >=f) ivs[nIv-1]->f=f;
		ivs[nIv-1]->t=t;
	}
	else{
		MapRange *iv=new MapRange(f,t);
		if(nIv >= capacity){
			capacity*=2;
			ivs=(MapRange**)realloc(ivs,capacity*sizeof(MapRange*));
		}
		ivs[nIv++]=iv;
	}
}


void IVSet::fin(){
	totLength=0;


	for(int i=0; i<nIv; i++){
		if(i > 0 && ivs[i]->f < ivs[i-1]->t) {
			ivs[i]->f = ivs[i-1]->t;
		}
		int l=ivs[i]->t - ivs[i]->f;
		totLength += l;
		ivs[i]->cumLength=totLength;
	}
	verb("Mapping:  lprofile=%i  totalIvLength=%i niv=%i\n",profileLength, totLength, nIv);
}
//===============================================================================


//===============================================================================
void IVSet::write(FILE*f){
	for(int i=0; i<nIv; i++){
		fprintf(f,"%i .. %i   %i   %i\n",ivs[i]->f,ivs[i]->t,ivs[i]->t-ivs[i]->f,ivs[i]->cumLength);
	}
}


int IVSet::randPos(){		// get a random position in the intervals
	int p=randInt(totLength);
	int i0=0, i1=nIv-1;
	int ivT=ivs[i1]->cumLength;
	int ivF=ivT-(ivs[i1]->t-ivs[i1]->f);
	if(ivF <= p && ivT >= p){
//		ivNo=i1;
		int pp=p-ivF; 		// position of rnd point in the interval
		pp=ivs[i1]->f+pp;	// position of rnd on the genome
		return pp;
	}


	while(1){				// binary search
		int i=(i0+i1)/2;
		ivT=ivs[i]->cumLength;
		ivF=ivT-(ivs[i]->t-ivs[i]->f);
		if(ivF <= p && ivT >= p){
//			ivNo=i;
			int pp=p-ivF; 		// position of rnd point in the interval
			pp=ivs[i]->f+pp;	// position of rnd on the genome
			return pp;
		}
		else if(p > ivT) i0=i+1;
		else if(p <=ivF) i1=i;
		if(i0==i1) break;
	}
	return 0;
}


void IVSet::print(int f, int t){
	for(int i=f; i<t && i<nIv; i++)
		printf("$iv$ %i..%i\n",ivs[i]->f,ivs[i]->t);
}





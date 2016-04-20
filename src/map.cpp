/*
 * map.cpp
 *
 *  Created on: 09.03.2013
 *      Author: 1
 */
#include "track_util.h"

MapPos::MapPos(){fg=0; pos=0; scaled=false;}
MapIv::MapIv(){beg.fg=true; end.fg=false;}

int MapPos::read(char* s){
	if(     *s=='B') fg=true;
	else if(*s=='E') fg=false;
	else return 1;
	pos=atoi(s+1);
	return 0;
}

int MapIv::read(char* s){
	char b[80];
	s=strtoupper(skipSpace(strcpy(b,s)));
	char *s1=s, *s2=strstr(s,"..");
	if(s2==0) return 1;
	*s2=0; s2+=2;
	if(beg.read(s1)) return 1;
	if(end.read(s2)) return 2;
	return 0;
}
void MapIv::scale(int pstep){beg.scale(pstep); end.scale(pstep);}

void MapPos::scale(int pstep){
	if(scaled) return;
	scaled=true;
	pos=(pos+pstep-1)/pstep;
}

void MapPos::print(){
	printf("%c:%i",fg?'b':'e',pos);
}

void MapPos::print(char *b){
	sprintf(b,"%c:%i",fg?'b':'e',pos);
}

void MapIv::print(){
	beg.print(); printf(".."); end.print();printf("\n");
}

char* MapIv::print(char *b){
	char bb[80], be[80];
	beg.print(bb); end.print(be);
	sprintf(b,"%s..%s",bb,be);
	return b;
}

//===============================================================================
void makeMapIv(bool cmpl, MapRange *mi);
void mapIntervals(bool cmpl, IVSet &ivs){
	for(int i=0; i<ivs.nIv; i++){
		makeMapIv(cmpl,ivs.ivs[i]);
	}
}
void mapIntervals(){
	mapIntervals(false,mapTrack.ivs);
//	deb("!!!! Map !!!!!");
//	for(int i=0; i<10; i++)	deb("%i..%i",mapTrack.ivs.ivs[i]->f,mapTrack.ivs.ivs[i]->t);
	if(mapTrack.hasCompl) mapIntervals(true,mapTrack.ivsC);
}

void fillMap(unsigned char *b, IVSet &ivs){
	zeroMem(b,profileLength);
	for(int i=0; i<ivs.nIv; i++){
		memset(b+ivs.ivs[i]->f, 1, ivs.ivs[i]->t-ivs.ivs[i]->f);
	}
}

void fillMap(){
	if(mapTrack.bytes==0){
		mapTrack.hasCompl=false;
		getMem(mapTrack.bytes,profileLength, "fill Map #1");
		memset(mapTrack.bytes, 1, profileLength);
		return;
	}
	fillMap(mapTrack.bytes, mapTrack.ivs);
	if(mapTrack.hasCompl)
		fillMap(mapTrack.cbytes, mapTrack.ivsC);
}
//======================= transform from-to to interval using boundaries
void makeMapIv(bool cmpl, MapRange *mi){
	int ff,tt;
	int f=mi->f, t=mi->t;
	if(!cmpl){
		if(miv.beg.fg) 	ff=f+miv.beg.pos;
		else 			ff=t+miv.beg.pos;
		if(miv.end.fg) 	tt=f+miv.end.pos;
		else 			tt=t+miv.end.pos;
	}
	else{
		if(miv.beg.fg) 	tt=t-miv.beg.pos;
		else 			tt=f-miv.beg.pos;
		if(miv.end.fg) 	ff=t-miv.end.pos;
		else 			ff=f-miv.end.pos;
	}
	if(ff<0) ff=0;
	mi->f=ff; mi->t=tt;
}

//===============================================================================
MapRange::MapRange(){f=t=cumLength=0;}
MapRange::MapRange(int fr, int to){ f=fr; t=to; cumLength=0;}

//===============================================================================
IVSet::IVSet(){
	capacity=1000; nIv=0; totLength=0; ivNo=0;
	errStatus="ivset init";
	getMem0(ivs,capacity, "ivset init #1");
	errStatus=0;
}

void IVSet::clear(){
	nIv=0; totLength=0;
	errStatus="ivset init";
	errStatus=0;
}
//===============================================================================
void IVSet::addIv(int f, int t){
//	f-=wProfSize; t+=wProfSize;
	if(f < 0) f = 0;    if(t > profileLength) t = profileLength;
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
		if(i > 0 && ivs[i]->f < ivs[i-1]->t) ivs[i]->f = ivs[i-1]->t;
		totLength += ivs[i]->t - ivs[i]->f;
		ivs[i]->cumLength=totLength;
	}
verb("Mapping:  lprofile=%i  totL=%i niv=%i\n",profileLength, totLength, nIv);
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
		ivNo=i1;
		int pp=p-ivF; 		// position of rnd point in the interval
		pp=ivs[i1]->f+pp;	// position of rnd on the genome
		return pp;
	}

	while(1){				// binary search
		int i=(i0+i1)/2;
		ivT=ivs[i]->cumLength;
		ivF=ivT-(ivs[i]->t-ivs[i]->f);
		if(ivF <= p && ivT >= p){
			ivNo=i;
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



/*
 * projector.cpp
 *
 *  Created on: 03 Jan 2017
 *      Author: Mironov
 */
#include "track_util.h"


char confTrackPath[4096];
char confProfPath[4096];
char confResPath[4096];
char confDir[1024];


bTrack *confBTrack;
bTrack *curTrack;
bTrack *projT;


FILE *prjLog=0;


void makeConfDir(char *confPath, char *path, const char*name, FILE *cfg){
	makeFileName(confPath,path,confDir);
	makeDir(confPath);
	fprintf(cfg,"%s=%s\n",name, confPath);
	strcat(confPath,"/");
}


void makeProj(char *fname){
	double xa=0,aa=0;
	char b[4096];
	if(projT->name) free(projT->name);
	verb("Projection: %s",fname);
	char *pp=profPath; profPath=confProfPath;	// set new path
	projT->initProfile(fname);
	profPath=pp;								// restore path


	for(int i=0; i<profileLength; i++){
		if(curTrack->isNA(i,false) || confBTrack->isNA(i,false)) continue;
		double a=confBTrack->getValue(i,false);
		double x=curTrack->getValue(i,false);
		xa+=a*x; aa+=a*a;
	}
	double z=xa/aa, e=0; int nn=0;
	double pMin=1.e+100, pMax=-1.e+100;
	for(int i=0; i<profileLength; i++){
		if(curTrack->isNA(i,false) || confBTrack->isNA(i,false)) continue;
		double a=confBTrack->getValue(i,false);
		double x=curTrack->getValue(i,false);


		double w=x-z*a;
		fProfile->set(i,w);
		e+=w; nn++;
		if(w < pMin) pMin=w;
		if(w > pMax) pMax=w;
	}


	verb("   min=%f  max=%f e=%f  n=%i z=%f\n",pMin, pMax, e/nn, nn, z);
	fprintf(prjLog,"<%s>\t%f\t%f\t%f\tz=%f\n",fname, pMin, pMax, e/nn,z);


	//=================================================== Write wig
	if(outPrjBGr){
		makeFileName(b,confTrackPath, fname, "bgraph");
		FILE *f=fopen(b,"wt");
		fprintf(f,"track type=wiggle_0 ");
		fprintf(f,"description=\"%s_projection_%s\" \n",curTrack->name,confFile);
		verb("Write profile %s...\n", curTrack->name);
		writeBedGr(f, fProfile);
		verb("\n");
		fclose(f);
	}


	projT->trackType=BED_GRAPH;
	projT->hasCompl=0;


	projT->finProfile();


	projT->writeProfilePrm(confProfPath);
	projT->writeByteProfile();
}


char currFname[4096];
void Projector(){
	confBTrack	=new bTrack();
	curTrack	=new bTrack();
	projT	=new bTrack();


	if(confFile==0) errorExit("Confounder not defined");
	if(fProfile==0) fProfile=new FloatArray();
	//================================= Make Directories
	char b[4096];
	sprintf(b,"%s.cfg",confFile);
	sprintf(confDir,"%s.proj",confFile);
	FILE * cfg=gopen(b,"wt");
	if(cfgFile) fprintf(cfg,"cfg=%s\n",cfgFile);
	makeConfDir(confTrackPath,trackPath, "trackPath", cfg);
	makeConfDir(confProfPath ,profPath , "profPath" , cfg);
	makeConfDir(confResPath  ,resPath  , "resPath"  , cfg);
	sprintf(b,"%s.bgraph",confFile);
	prepare(b);
//	int cxcx=confBTrack->openTrack(b);
	prjLog=fopen("projections","w");
	fprintf(prjLog,"confounder=<%s>\n",confFile);
	fprintf(prjLog,"track\tminVal\tmaxVal\te\tproj_coef\n");
	for(int i=0; i<nfiles; i++){
		strcpy(currFname, files[i].fname);
		prepare(files[i].fname);
		curTrack->openTrack(files[i].fname);
		makeProj(files[i].fname);
		curTrack->clear();
		projT->clear();
	}
	fclose(prjLog);
	del(fProfile); fProfile=0;
}





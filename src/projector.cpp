/*
 * projector.cpp
 *
 *  Created on: 03 Jan 2017
 *      Author: Mironov
 */
#include "track_util.h"


char confTrackPath[TBS];
char confProfPath[TBS];
char confResPath[TBS];
char confDir[TBS];


bTrack *confBTrack;
bTrack *curTrack;
bTrack *projT;


FILE *prjLog=0;


void makeConfDir(char *confPath, const char *path, char * tPath, char *tName, const char*prmName, FILE *cfg){
	sprintf(confPath, "%s%s%s.proj/",path,tPath,tName);	//========== additional path in the track
	makeDir(confPath);
}


void makeProj(){
	double xa=0,aa=0;
	char b[TBS];
//	if(projT->name) free(projT->name);
	verb("Projection: %s",curTrack->name);
	char *pp=profPath; profPath=confProfPath;	// set new path
	projT->initProfile(curTrack->name);
	profPath=pp;								// restore path

//==========	calculate scalar productions
	for(int i=0; i<profileLength; i++){
		if(curTrack->isNA(i,false) || confBTrack->isNA(i,false)) continue;
		double a=confBTrack->getValue(i,false);
		double x=curTrack->getValue(i,false);
		xa+=a*x; aa+=a*a;
	}
	double z=xa/aa, e=0; int nn=0;
	double pMin=1.e+100, pMax=-1.e+100;
//=========== Make projections
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
	fprintf(prjLog,"<%s>\t%f\t%f\t%f\tz=%f\n",curTrack->name, pMin, pMax, e/nn,z);


	//=================================================== Write track
	if(outPrjBGr){
		curTrack->makePath(confTrackPath);
		curTrack->makeFname(b,confTrackPath,"bgraph");
		FILE *f=xopen(b,"wt");
		fprintf(f,"track type=wiggle_0 ");
		fprintf(f,"description=\"%s_projection_%s\" \n",curTrack->name,confounder);
		verb("Write profile %s...\n", curTrack->name);
		writeBedGr(f, fProfile);
		verb("\n");
		fclose(f);
	}


	projT->trackType=BED_GRAPH;
	projT->hasCompl=0;

	projT->finProfile();

	projT->writeProfilePrm();
	projT->writeByteProfile();
}

void copyCfg(FILE *cfg){
	char b[TBS];
	FILE *old=fopen(cfgFile,"r");
	if(old==0) return;
	for(; fgets(b,sizeof(b),old);){
		char bb[TBS];
		char *s=trim(strcpy(bb,b));
		char *ss=strchr(s,'=');
		if(ss){
			*ss=0; ss=trim(s);
			if(keyCmp(s,"confounder") == 0) {fprintf(cfg,"#%s\n",b); continue;		}
			if(keyCmp(s,"trackPath") == 0) {fprintf(cfg,"trackPath=%s\n",confTrackPath); continue;}
			if(keyCmp(s,"profPath" ) == 0) {fprintf(cfg,"profPath=%s\n" ,confProfPath  );continue;}
			if(keyCmp(s,"resfPath" ) == 0) {fprintf(cfg,"profPath=%s\n" ,confProfPath  );continue;}
		}
	fprintf(cfg,"%s",b);
	}
	fprintf(cfg,"\n\n");
	fclose(old);
}

char currFname[TBS];
void Projector(){
	if(fileExists(defaultConfig))	//referense to existing profile file
		cfgFile=strdup(defaultConfig);

	confBTrack	=new bTrack(confounder);
	curTrack	=new bTrack();
	projT	=new bTrack();

	if(fProfile==0) fProfile=new FloatArray();
	//================================= Make Directories
	char b[TBS];
	confBTrack->makeFname(b, (char*)("./"),"cfg");
	FILE * cfg=xopen(b,"wt");
	snprintf(confDir,sizeof(confDir), "%s.proj",confBTrack->name);

	makeConfDir(confTrackPath, trackPath, confBTrack->path, confBTrack->name , "trackPath", cfg);
	makeConfDir(confProfPath , profPath, confBTrack->path , confBTrack->name, "profPath" , cfg);
	makeConfDir(confResPath  ,resPath, confBTrack->path  , confBTrack->name, "resPath"  , cfg);

	if(fileExists(cfgFile)) copyCfg(cfg);

	profPath =confProfPath;
	prjLog=fopen("projections","w");
	fprintf(prjLog,"confounder=<%s>\n",confounder);
	fprintf(prjLog,"track\tminVal\tmaxVal\te\tproj_coef\n");


	for(int i=0; i<nfiles; i++){
		strcpy(currFname, files[i].fname);
		prepare(files[i].fname);
		curTrack->openTrack(files[i].fname);
		makeProj();
		curTrack->clear();
		projT->clear();
	}
	fclose(prjLog);
	del(fProfile); fProfile=0;
}





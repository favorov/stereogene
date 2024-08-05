/*
 * covar_mtx.cpp
 *
 *  Created on: 24 May 2017
 *      Author: Mironov
 */
#include "track_util.h"


bTrack **tracks;
CovarMtx cMtx;
VectorX *eVector=0;
double   eValue=0;
char* confFile=0;


void Confounder(){
	bTrack *confdr=new bTrack();
	if(confFile==0) confFile=strdup("confounder");
	confdr->trackType=BED_GRAPH;
	confdr->hasCompl=0;
	confdr->name=confFile;
	confdr->initProfile();

	fProfile->init(0);

	FILE *f;
	char b[TBS];
	verb("\n");
	verb("Make Confounder...\n");
	for(int j=0; j<nfiles; j++){
		for(int i=0; i<profileLength; i++){
			if(i%1000000 ==0) verb("%5.1f%%\r",1.*i/profileLength*100);
			double w=tracks[j]->getValue(i,0);
			if(w==NA) w=0;
//			w=(w-tracks[j]->nativeAv)/tracks[j]->nativeSd; //normalize track by Z-score
			w=(w)/tracks[j]->nativeAv;	//normalize track by average
			double x=fProfile->get(i);
			x+=eVector->v[j]*w;
			fProfile->set(i,x);
		}
	}
	verb("\n");
	//=================================================== Norm
	confdr->finProfile();
	//=================================================== Write bedgraph
//	makeFileName(b,sizeof(b), trackPath, confFile, BGR_EXT);
	makeFileName(b, trackPath, confFile, BGR_EXT);
	f=xopen(b,"wt");
	fprintf(f,"track type=bedGraph name=\"%s\" ", confFile);
	fprintf(f,"description=\"confounder\" \n");
	for(int i=0; i<nfiles; i++){
		fprintf(f,"#%s\n",files[i].fname);
	}
	verb("Write confounder profile...\n");

	writeBedGr(f,fProfile);
	verb("\n");
	fclose(f);

	confdr->writeProfilePrm();
	confdr->writeByteProfile();
	del(fProfile); fProfile=0;
	del(confdr);
}


//===================================================================
//===================================================================


CovarMtx::CovarMtx(int nn){
	init(nn);
}
void CovarMtx::init(int nn){
	Matrix::init(nn);
	int n2=nn*nn;
	getMem(cov  ,n2 ,"covar1");	 zeroMem(cov,  n2);
	getMem(mean ,n  ,"covar2");  zeroMem(mean,n);
	getMem(sd   ,n  ,"covar2");  zeroMem(sd,n);
	getMem(count,n  ,"covar3");  zeroMem(count,n2);
}


CovarMtx::~CovarMtx(){
	if(cov  ) xfree(cov  , "free covar 1");
	if(mean ) xfree(mean , "free covar 2");
	if(sd   ) xfree(sd   , "free covar 3");
	if(count) xfree(count, "free covar 4");
}


//===================================================================================
//================== Calculate the covariations for pair of the tracks ==============
void CovarMtx::addCov(int itrack, int jtrack){
	Track *tr1=tracks[itrack];
	Track *tr2=tracks[jtrack];
	double c=0, e=0, d=0, n=0;


	for(int i=0; i<profileLength; i++){
		if(tr1->isNA(i,0) || tr1->isNA(i,0) ) continue;
		double x1=tr1->getValue(i,0);
		double x2=tr2->getValue(i,0);

		if(itrack==0) {e+=x2; d+=x2*x2;}
		c+=x1*x2;
		n++;
	}
	int idx=getIndex(itrack,jtrack);

	if(itrack==0){
		mean [jtrack]=e;
		sd   [jtrack]=d;
		count[jtrack]=n;
	}
	cov  [idx   ]=c;

}


//===================================================================================
double CovarMtx::calc(int itrack, int jtrack){

	addCov(itrack,jtrack);

	int idx=getIndex(itrack,jtrack);
	double c=cov[idx], e1=mean[itrack], e2=mean[jtrack], d1=sd[itrack], d2=sd[jtrack];
	int n1=count[itrack], n2=count[jtrack];

	e1/=n1; e2/=n2;
	d1=sqrt ( (d1-n1*e1*e1)/(n1-1));
	d2=sqrt ( (d2-n2*e2*e2)/(n2-1));

	double n12=sqrt(double(n1)*double(n2));
	c=(c-e1*e2*n12)/n12;
	c=c/(d1*d2);

	values[getIndex(itrack,jtrack)]=
			values[getIndex(jtrack,itrack)]=c;
	verb(" cov[%i,%i]=%f\n",itrack,jtrack,c);
	return c;
}


//===================================================================================
void CovarMtx::print(FILE *f){
	for(int i=0; i<nfiles; i++){
		fprintf(f,"\t%s",tracks[i]->name);
	}
	fprintf(f,"\n");
	for(int i=0; i<nfiles; i++) {
		fprintf(f,"%s",tracks[i]->name);
		for(int j=0; j<nfiles; j++){
			fprintf(f,"\t%.4f",get(i,j));
		}
	fprintf(f,"\n");
	}
	fflush(f);
}


//===========================================================


void calcCovar(){
	CovarMtx *cMtx=new CovarMtx(nfiles);

	for(int i=0; i<nfiles; i++){
		for(int j=i; j<nfiles; j++){
			verb("Covariations %i~%i: <%s>~<%s> ",i,j,tracks[i]->name,tracks[j]->name);
			cMtx->calc(i,j);
		}
	}
	eVector=new VectorX(nfiles);
	Matrix *mtx=new Matrix(cMtx);
	eValue = mtx->eigen(eVector);
	char b[TBS];
	sprintf(b,"%s.covar",confFile);
	FILE *f=xopen(b,"wt");
	cMtx->print(f);

	fputs("eigenVector =",f);
	eVector->print(f);
	fputs("\n",f);
	fprintf(f,"eigenValue=%.5f; \n",eValue);
	fclose(f);
	del(cMtx);
	del(mtx);
}


//================================================================================
//================================================================================
//================================================================================
//================================================================================
//================================================================================


void Covariator(){
	getMem(tracks, nfiles,"");
	for(int i=0; i<nfiles; i++){
		tracks[i]=new bTrack();
		verb("read %s\n",files[i].fname);
		tracks[i]->openTrack(files[i].fname);
	}
	if(confFile==0) confFile=strdup("confounder");
	calcCovar();
	Confounder();


	for(int i=0; i<nfiles; i++) del(tracks[i]);
	xfree(tracks,"covariator: free tracks");
	exit(0);
}





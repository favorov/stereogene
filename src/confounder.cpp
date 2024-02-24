/*
 * covar_mtx.cpp
 *
 *  Created on: 24 May 2017
 *      Author: Mironov
 */
#include "track_util.h"


bTrack **tracks;
CovarMtx cMtx;
Matrix *eVectors=0;
double *eValues=0;


//=============== SG-covariator  =================
//================================================
//int qmin=98,qmax=101;
VectorX::VectorX(int nn){init(nn);}
VectorX::VectorX(){init(nfiles);}
void VectorX::init(int nn){	n=nn; getMem(v,n,""); zeroMem(v,n);}
//==========================================================================================
void VectorX::get(int pos){
	get(pos,v);
}
//==========================================================================================
void VectorX::get(int pos, double *val){
	for(int i=0; i<nfiles; i++){
		if(tracks[i]->isNA(pos,0))val[i]=NA;
		else val[i]=tracks[i]->getValue(pos,0);
	}
}


double VectorX::scalar(VectorX *x){
	return scalar(*x);
}
double VectorX::scalar(VectorX &x){
	double w=0;
	int nn=0;
	for(int i=0; i<n; i++) {
		if(v[i]==NA || x.v[i]==NA) continue;
		w+=v[i]*x.v[i]; nn++;
	}
	if(nn==0) return NA;
	w=w/nn*n;
	return w;
}


//==========================================================================================
int VectorX::chk(int pos){
	double b[n];
	get(pos, b);
	int nz=0;
	for(int i=0; i<n; i++) if(b[i]<=1) nz++;
	return nz;
}
//==========================================================================================
int VectorX::chk(){
	int nz=0;
	for(int i=0; i<n; i++) if(v[i]<=1) nz++;
	return nz;
}
//==========================================================================================
void VectorX::printH(FILE *f){
	for(int i=0; i<nfiles; i++){
		fprintf(f,"%s", tracks[i]->name);
		if(i<n-1) fputc('\t',f);
		else	  fputc('\n',f);
	}
}


//==========================================================================================
void VectorX::print(FILE *f){
	for(int i=0; i<n; i++){
		fprintf(f,"%.3f",v[i]);
		if(i<n-1) fputc('\t',f);
		else	  fputc('\n',f);
	}
}
//==========================================================================================
//==========================================================================================


void Confounder(){
	bTrack *confdr=new bTrack();
	if(confFile==0) confFile=strdup("confounder");
	confdr->trackType=BED_GRAPH;
	confdr->hasCompl=0;
	confdr->name=confFile;
	confdr->initProfile();
	VectorX v=VectorX();
	VectorX cnf=VectorX();
	FILE *f;
	char b[4096];
	for(int i=0; i<nfiles; i++){
		cnf.v[i]=eVectors->get(i,0);
	}
	verb("\n");
	if(verbose) {verb("confounder: "); cnf.print(stdout);}


	double min=1.e+100, max=-1.e+100, e=0;
	int nn=0;
	int na0=0;
	verb("Make Confounder...\n");
	for(int i=0; i<profileLength; i++){
		if(i%1000000 ==0) verb("%5.1f%%\r",1.*i/profileLength*100);
		v.get(i);
		if(v.v[0]==NA) na0++;
		double w=v.scalar(cnf);
		if(w==NA){
			fProfile->set(i,w);
			continue;
		}
		fProfile->set(i,w);	//finProfile will divide by bisize
		e+=w; nn++;
		if(w<min) min=w;
		if(w>max) max=w;
	}
	verb("\n");
	e/=nn;
	if(e < 0){
		min=1.e+100; max=-1.e+100; e=0;
		for(int i=0; i<profileLength; i++){
			double w=fProfile->get(i);
			if(w==NA) continue;
			fProfile->set(i,(w=-w));
			e+=w;
			if(w<min) min=w;
			if(w>max) max=w;
		}
	}
	//================================================ Norm
	confdr->finProfile();
	//=================================================== Write wig
	makeFileName(b,trackPath, confFile, BGR_EXT);
	f=gopen(b,"wt");
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
	getMem(cov,	 n2,"covar1");	 zeroMem(cov,  n2);
	getMem(meani,n2,"covar2");   zeroMem(meani,n2);
	getMem(meanj,n2,"covar2");   zeroMem(meanj,n2);
	getMem(count,n2,"covar3"); 	 zeroMem(count,n2);
}


CovarMtx::~CovarMtx(){
	if(cov)   xfree(cov, "free covar 1");
	if(meani) xfree(meani, "free covar 2");
	if(meanj) xfree(meanj, "free covar 3");
	if(count) xfree(count, "free covar 3");
}


//===================================================================================
void CovarMtx::addCov(int itrack, int jtrack, int f, int t){
	Track *tr1=tracks[itrack];
	Track *tr2=tracks[jtrack];
	double c=0, e1=0, e2=0, n=0;


	for(int i=f; i<t; i++){
		if(tr1->isNA(i,0) || tr1->isNA(i,0) ) continue;
		double x1=tr1->getValue(i,0);
		double x2=tr2->getValue(i,0);


		if(itrack==0) e2+=x2;
		if(jtrack==0) e1+=x1;
		c+=x1*x2;
		n++;
	}
	int idx=getIndex(itrack,jtrack);
	meani[idx]+=e1;
	meanj[idx]+=e2;
	cov  [idx]+=c;
	count[idx]+=n;
}




//===================================================================================
double CovarMtx::calc(int itrack, int jtrack){
	int l=profileLength;
	addCov(itrack,jtrack, 0, l);
	int idx=getIndex(itrack,jtrack);
	double c=cov[idx], e1=meani[idx], e2=meanj[idx];
	int n=count[idx];
	e1/=n; e2/=n;
	c-=e1*e2*n;
	c/=n;
	values[getIndex(itrack,jtrack)]=
			values[getIndex(jtrack,itrack)]=c;
	verb(" cov=%f\n",c);
	return c;
}


//===================================================================================
void CovarMtx::print(FILE *f){
	for(int i=0; i<nfiles; i++) fprintf(f,"\t%s",tracks[i]->name);
	fprintf(f,"\n");
	for(int i=0; i<nfiles; i++) {
		fprintf(f,"%10s",tracks[i]->name);
		for(int j=0; j<nfiles; j++){
			fprintf(f,"\t%.4f",get(i,j));
		}
	fprintf(f,"\n");
	}
}


//===========================================================


void calcCovar(){
	int nIter=100;
	cMtx.init(nfiles);


	for(int i=0; i<nfiles; i++){
		for(int j=i; j<nfiles; j++){
			verb("Covariations %i~%i: <%s>~<%s> ",i,j,tracks[i]->name,tracks[j]->name);
			cMtx.calc(i,j);
		}
	}
	char b[2048];
	sprintf(b,"%s.cvr",confFile);
	FILE *f=xopen(b,"wt");
	cMtx.print(f);
	Matrix *x = new Matrix(&cMtx);
	getMem(eValues,nfiles+10,"eigenVal"); zeroMem(eValues,nfiles);
	eVectors=eigenVectors(x,eValues,nIter,1.E-3);
	fputs("eigenVectors\n",f);
	eVectors->printMtx(f);
	fputs("eigenValues\n",f);
	for(int i=0; i<nfiles; i++)
		fprintf(f,"%.5f; ",eValues[i]);
	fprintf(f,"\n");
	fclose(f);
	del(x);
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
	calcCovar();
	Confounder();


	for(int i=0; i<nfiles; i++) del(tracks[i]);
	xfree(eValues,"covariator: free eValues");
	xfree(tracks,"covariator: free tracks");
	exit(0);
}





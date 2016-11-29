/*
 * bTrack.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: mironov
 */
#include "track_util.h"


void removeNA(unsigned char *s, int l){
	if(NAFlag) return;
	for(int i=0; i<l; i++){
		if(s[i]==0) s[i]=1;
		if(inpThreshold){
			if((s[i]-1)*100./250 > inpThreshold) s[i]=250;
			else s[i]=1;
		}
	}
}

void failParam(const char* s){
	verb("parameter \'%s\' failed\n",s);
}

int chkVersion(char *ver){
	char ver1[1024], ver0[1024];
	return strcmp(getMajorVer(version,ver0), getMajorVer(ver,ver1));
}

bool bTrack::check(const char *fname){
	if(clearProfile){
		verb("forced profile recalculation\n");
		return false;
	}
	char prmFile[4096], binFile[4096], b[4096];
	makeFileName(prmFile,profPath, name, PRM_EXT);
	makeFileName(binFile,profPath, name, BPROF_EXT);
//======================================================= check existing files
	if(! fileExists(binFile)) {
		verb("file < %s > does not exist\n",binFile); return false;
	}
	if(! fileExists(prmFile)) {
		verb("file < %s > does not exist\n",prmFile); return false;
	}
	//================= Check modification time
	makeFileName(b,trackPath, name);
	unsigned long tPrm  =getFileTime(prmFile);
	unsigned long tTrack=getFileTime(b);

	if(tPrm < tTrack){
		verb("Old profile prm_time=%li  track_time=%li\n",tPrm,tTrack); return false;
	}

//======================================================= check parameters
	FILE* f=xopen(prmFile,"rt");
	trackType=getTrackType((const char*)name);
	bool fg=true;
	char bver[80]="";
	for(;fgets(b,sizeof(b),f)!=0;){
		char *s1=strtok(b,"=");
		char *s2=strtok(0,"=\r\n");
		if(b[0]=='#') continue;

		else if(strcmp(s1,"version")==0){
			strcpy(bver,s2);
		}
		if(strcmp(s1,"input"    )==0){
			if(trackType == MODEL_TRACK){
				model.readMap(name);
				if(strcmp(s2,model.definition)) {failParam("input");  fg=false; break;}
			}
			else if(strcmp(s2,name)) {
				failParam("input");  fg=false; break;}}
		else if(strcmp(s1,"step")==0){
			if(atoi(s2)!=binSize) {failParam("step");  fg=false; break;}
		}
		else if(strcmp(s1,"bpType")==0){
			if(atoi(s2)!=bpType) {failParam("bpType");  fg=false; break;}
		}
		else if(strcmp(s1,"intervFlag")==0){
			if(atoi(s2)!=intervFlag0) {failParam("intervFlag");  return false;}
		}
		else if(strcmp(s1,"NA")==0){
			if(atoi(s2)!=NAFlag) {failParam("intervFlag");  return false;}
		}
	}
	fclose(f);
	if(chkVersion(bver)){
		failParam("program version"); return false;
	}
	return fg;
}


//=============================================================================
void bTrack::read(const char *fname){

	for(;*fname==DERIV; fname++) {
		deriv++;
	}
//todo  If this is BED and ivFlag defined then define extended filename
	name=strdup(fname);
	char prmFile[4096], binFile[4096], b[4096];
	makeFileName(prmFile,profPath, name, PRM_EXT);
	makeFileName(binFile,profPath, name, BPROF_EXT);
	if(!check(fname)) makeBinTrack();

	//============================================== Read params
	errStatus="read bytes";
	FILE *f=xopen(prmFile,"rt");
	for(;fgets(b,sizeof(b),f)!=0;){
		char *s1=strtok(b,"=");
		char *s2=strtok(0,"=\r\n");
		if(     strcmp(s1,"min"    )==0) minP=atof(s2);
		else if(strcmp(s1,"max"    )==0) maxP=atof(s2);
		else if(strcmp(s1,"average")==0) av0 =atof(s2);
		else if(strcmp(s1,"stdDev" )==0) sd0 =atof(s2);
		else if(strcmp(s1,"strand" )==0) hasCompl=atoi(s2);
		else if(strcmp(s1,"lscale" )==0) lScale  =atoi(s2);
		else if(strcmp(s1,"scaleFactor" )==0) scaleFactor=atof(s2);
		else if(strcmp(s1,"intervFlag"  )==0) intervFlag =atoi(s2);
	}
	if(sd==NAN) sd=1;
	fclose(f);
	//============================================== Read bytes
	f=xopen(binFile,"rb"); if(f==0) return;

	//============================================= get file length
	fseek(f, 0L, SEEK_END); lProf = (int) ftell(f);
	fseek(f, 0L, SEEK_SET);
	if(hasCompl) lProf/=2;
	getMem0(bytes,(lProf+10), "Read bTrack");

	fread(bytes,1,lProf,f); removeNA(bytes,lProf);
	if(hasCompl){
		getMem0(cbytes,(lProf+10), "bTrack read #1");
		fread(cbytes,1,lProf,f); removeNA(cbytes,lProf);
	}
	fclose(f);

	delta=minP;
	bScale=(maxP-minP)/250.;

	av=sd=nn=0;

	getMem0(profWindow,(profWithFlanksLength+10), "bTrack read #2");
	errStatus=0;
//	if(lAuto) trackAutoCorrelation();
}

//===============================================================================
double bTrack::addStatistics(){
	double avv=0, sdd=0;
	for(int i=0; i<profWithFlanksLength; i++){
		double x=profWindow[i];
		avv+=x; sdd+=x*x;  nn++;
	}
	av+=avv; sd+=sdd;
	return avv/profWithFlanksLength;
}
void bTrack::finStatistics(){
	av/=nn; sd=sd-av*av*nn; sd/=nn-1; sd=sqrt(sd);
}

//========================================================================
void multByteProf(unsigned char *bytes,unsigned char *p1,unsigned char *p2){
	if(p1==0 && p2==0) return;
	for(int i=0; i<profileLength; i++){
		int x=0;
		if(p1) x|=p1[i];
		if(p2) x|=p2[i];
		bytes[i]=bytes[i]*x;
	}
}
void bTrack::makeIntervals(){
	if(hasCompl) {
		multByteProf(bytes,mapTrack.bytes,0);
		makeIntervals(bytes, &ivs);
		if(mapTrack.hasCompl){
			multByteProf (cbytes,mapTrack.cbytes,0);
			makeIntervals(cbytes,&ivsC);
		}
		else{
			multByteProf (cbytes, mapTrack.bytes,0);
			makeIntervals(cbytes, &ivsC);
		}
	}
	else if(mapTrack.hasCompl){
		multByteProf (bytes, mapTrack.bytes, mapTrack.cbytes);
		multByteProf (bytes, mapTrack.bytes, 0);
		makeIntervals(bytes, &ivs);
	}
	else{
		multByteProf(bytes,mapTrack.bytes,0);
		makeIntervals(bytes, &ivs);
	}

	if(ivs.nIv+ivsC.nIv == 0) errorExit("no nonZero windows");
	ivs.fin();
	if(ivsC.nIv) ivsC.fin();

	int l0=lProf; if(hasCompl) l0*=2;
	int l1=ivs.totLength; if(hasCompl) l1+=ivsC.totLength;
	av0=av0*l0/l1;
	sd0=sd0*sqrt(l0/l1);

}
//========================================================================

int isZero(unsigned char *bytes, int i){
	if(bytes[i] <2 || bytes[i] <= threshold) return 1;
	return 0;
}
void bTrack::makeMapIntervals(){
	makeMapIntervals(bytes, &ivs);
	if(hasCompl) makeIntervals(cbytes,&ivsC);
//	deb("!======= Map =======");
//	for(int i=0; i<10; i++)	deb("%i..%i",ivs.ivs[i]->f,ivs.ivs[i]->t);
}

//========================================================================
void bTrack::makeMapIntervals(unsigned char *bytes, IVSet *iv){
	int f=0;
	bool space=true;
	for(int i=0; i<lProf; i++){
		if(isZero(bytes,i)){
			if(!space){
				iv->addIv(f,i);
			}
			space=true;
			}
		else{
			if(space) f=i;
			space=false;
		}
	}
	if(!space) iv->addIv(f,lProf);
}


void bTrack::makeIntervals(unsigned char *bytes, IVSet *iv){
	int f=0;
	bool space=true;

	int nZ=0;
	int wp=wProfSize;
	for(int i=0; i<wp; i++)
		if(isZero(bytes,i)) nZ++;
	for(int i=wp; i<lProf; i++){
		if(nZ >= wp-1){
			if(!space) iv->addIv(f,i-wp);
			space=true;
			}
		else{
			if(space) f=i-wp;
			space=false;
		}

		if( isZero(bytes,i          )) nZ++;
		if( isZero(bytes,i-wp)) nZ--;
	}
	if(!space) iv->addIv(f,lProf);
}

//========================================================================
int bTrack::getRnd(bool cmpl){
	int pos;
	pos=cmpl ? ivsC.randPos() : ivs.randPos();         // get random position in the interval
	return pos;
}


//========================================================================
void bTrack::clear(){
	ivs.clear();
}

//========================================================================
int  bTrack::countNA(int pos, bool cmpl){
	int c=pos+wProfSize >= lProf ? wProfSize+pos-lProf-1 : 0;
	for(int i=0; i< wProfSize && i+pos < lProf; i++) {
		if(complFg==IGNORE_STRAND){		// ignore strand
			int x=1;
			if(bytes[pos+i]!=0) x=0;
			if(hasCompl && cbytes[pos+i]==0) x=0;
			c+=x;
		}
		else if(!cmpl || !hasCompl) 		// colinear or profile do not know the orient
			{if(bytes[pos+i]==0) 	c++;}
		else if(cmpl)						// comlement
			{if(cbytes[pos+i]==0) 	c++;}
	}
	return c;
}

//========================================================================
int  bTrack::countZero(int pos, bool cmpl){
	int c=pos+wProfSize >= lProf ? wProfSize+pos-lProf-1 : 0;
	int thr= threshold;
	for(int i=0; i< wProfSize && i+pos < lProf; i++) {
		if(complFg==IGNORE_STRAND){		// ignore strand
			int x=1;
			if(bytes[pos+i] > thr) x=0;
			if(hasCompl && cbytes[pos+i] > thr) x=0;
				c+=x;
		}
		else if(cmpl || !hasCompl) 		// colinear or profile do not know the orient
			{if(bytes[pos+i] <= thr) 	c++;}
		else if(cmpl){					// comlement
			{if(cbytes[pos+i] <= thr) 	c++;}
		}
	}
	return c;
}

bTrack::bTrack(const char* fname){read(fname);}

bTrack::bTrack(){
	cbytes=bytes=0;name=0;profWindow=0;
	avWindow=sdWindow=0;// mean and stdDev in current window

	lProf=lScale=0;
	av=sd=minP=maxP=bScale=delta=nn=0;
	av0=sd0=0;
	strandFg=1; deriv=0;
	trackType=intervFlag=0;
	projCoeff=0;
	hasCompl=false;
	scaleFactor=0.2;
}
//======================================= decode binary value to real value
double bTrack::getVal(unsigned char b){
	if(b==0 && NAFlag && trackType==WIG_TRACK) return rGauss()*sd*noiseLevel;
	if(b < threshold) return 0;
	if(b==1) return 0;
	double x=bScale*(b-1)+delta;
	if(lScale) x=(exp(x)-1)/scaleFactor;
	return x;
}

//======================================= decode binary value at certain position
double bTrack::getValue(int pos, int cmpl){
	if(cmpl){
		if(hasCompl) return getVal(cbytes[pos]);
		else return getVal(bytes[pos]);
	}
	return getVal(bytes[pos]);
}

//====================================== calculate projection to orthogonal subspace
void bTrack::ortProject(){
	projCoeff=0;
	if(pcorProfile==0) return;
	double xy=0, xx=0;
	for(int i=0; i<lProf; i++){
		double x,y;
		x=projTrack.getValue(i,false);
		y=getValue(i,false);
		xx+=x*x; xy+=x*y;
		if(hasCompl){
			x=projTrack.getValue(i,true);
			y=getValue(i,true);
			xx+=x*x; xy+=x*y;
		}
	}
	projCoeff=xy/xx;
}

//================================================= Get projected value
double bTrack::getProjValue(int pos, bool cmpl){
	if(pos >= lProf) return 0;
	double x=getValue(pos,cmpl);
	if(projCoeff) x-=projCoeff*projTrack.getValue(pos,cmpl);
	return x;
}

//================================================= decode the values to an array
double * bTrack::getProfile(int pos, bool cmpl){ //====== pos - profile position; cmpl=true <=> +strand
	double *a=profWindow;
	//======================================================= fill window
	for(int i=0; i < wProfSize; i++, a++){
		double x=0,y=0;
		*a=0;
		if(complFg==IGNORE_STRAND){				// ignore strand
			x = getProjValue(pos+i+LFlankProfSize,false);		// profile
			if(hasCompl)
				y = getProjValue(pos+i+LFlankProfSize,true);				// + complement profile
		}
		else if(cmpl)                          // colinear or profile do not know the orient
			y=getProjValue(pos+i+LFlankProfSize,hasCompl);
		else
			x=getProjValue(pos+i+LFlankProfSize,false);
		*a=x+y;
		avWindow+=*a; sdWindow+=(*a) *(*a);
	}
	//======================================= Window Statistics
	avWindow/=wProfSize;
	sdWindow=sdWindow-avWindow*avWindow/wProfSize;
	sdWindow/=wProfSize-1; sdWindow=sqrt(sdWindow);
	//======================================== fill flanks
	int x0=wProfSize+LFlankProfSize;
	int x1=x0+LFlankProfSize+RFlankProfSize;
	for(int x=x0; x<x1; x++){
		double xq=rGauss(avWindow,sdWindow);
		profWindow[x%profWithFlanksLength]=xq;
	}
	return profWindow;
}
//======================================================
double *autoCorr;
double *dat;
int lProfAuto;

double arrayMax(double *d, int n){
	double mre=-1.e+200;
		for(int i=0; i<n; i++) mre=max(mre,d[i]);
	return mre;
}

////========================================== Generate autocorrelation function for an array
//void autoCorrelation(double *p, int from, int to){
//	Fourier ff(lProfAuto);
//	zeroMem(dat,lProfAuto);
//
//	for(int i=from; i<to; i++){
//		int k=(i-from)/corrScale;
//		if(p[k]!=NA && p[k]!=0) {
//			dat[k]+=p[k];
//		}
//	}
//
//	if(arrayMax(dat,lProfAuto) <=0) return;
//
//	ff.calc0(dat,0);
//	for(int i=0; i<lProfAuto; i++){
//		dat[i]=ff.re[i]*ff.re[i]+ff.im[i]*ff.im[i];
//	}
//
//	ff.calc0(dat,0);
//	double mre=arrayMax(ff.re, lProfAuto);
//	if(mre==0) return;
//
//	for(int i=0; i<lProfAuto; i++) autoCorr[i]+=ff.re[i]/mre;
//}
//
////========================================= Generate autocorrelation function for track
//void bTrack::trackAutoCorrelation(){
//	verb("Autocorrelation...\n");
//	errStatus="Autocorrelation";
//	corrScale=1;
//
//	lProfAuto=lAuto/(binSize*corrScale);
//
//	getMem(xDat,lProfAuto, "Track auto #1");
//	getMem(yDat,lProfAuto, "Track auto #2");
//
//	getMem0(autoCorr,lProfAuto, "Track auto #3"); zeroMem(autoCorr,lProfAuto);
//	getMem0(dat,lProfAuto, "Track auto #4");
//
//	int from=0, to=lProfAuto;
//	for(; to< profileLength; from+=lProfAuto, to+=lProfAuto){
//		verb("%i      \r",from);
//
//		if(readProfileToArray(xDat,corrScale,from,to,false)==0) return;
//		autoCorrelation(xDat,from,to);
//		if(profilec!=0){
//			if(readProfileToArray(yDat,corrScale,from,to,false)==0) return;
//			autoCorrelation(yDat,from,to);
//		}
//	}
//	verb("\n");
//
//	char b[1024];
//	//============================================= Write autocorrelation
//	double ma=-1.e+200;
//	for(int i=0; i<lProfAuto; i++) ma=max(autoCorr[i],ma);
//	for(int i=0; i<lProfAuto; i++) autoCorr[i]/=ma;
//	makeFileName(b,resPath,name,"ac");
//	FILE *fil=xopen(b,"wt");
//	fprintf(fil,"l\tac\n");
//	for(int i=lProfAuto/2; i>0; i--){
//		if(abs(autoCorr[i])<0.01) continue;
//		int k=binSize*(i-1);
//		fprintf(fil,"%i\t%.5f\n",-k,autoCorr[i]);
//	}
//	for(int i=1; i<lProfAuto/2; i++){
//		if(abs(autoCorr[i])<0.01) continue;
//		int k=binSize*(i-1);
//		fprintf(fil,"%i\t%.5f\n",k,autoCorr[i]);
//	}
//	fclose(fil);
//	xfree(xDat,"xDat");
//	xfree(yDat,"yDat");
//	xfree(autoCorr,"autoCorr");
//	xfree(dat,"dat");
//	verb("OK\n");
//	errStatus=0;
//}
//
//

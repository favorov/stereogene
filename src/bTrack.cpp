/*
 * bTrack.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: mironov
 */
#include "track_util.h"


void removeNA(BINVAL *s, int l){
	if(NAFlag) return;
	for(int i=0; i<l; i++){
		if(s[i]==NA) s[i]=0;
		if(inpThreshold){
			if((s[i])*100./250 > inpThreshold) s[i]=250;
			else s[i]=0;
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
	makeFileName(prmFile,profPath, fname, PRM_EXT);
	makeFileName(binFile,profPath, fname, BPROF_EXT);
//======================================================= check existing files
	if(! fileExists(binFile)) {
		verb("file < %s > does not exist. Preparing binary track\n",binFile); return false;
	}
	if(! fileExists(prmFile)) {
		verb("file < %s > does not exist. Preparing binary track\n",prmFile); return false;
	}
	//================= Check modification time
	makeFileName(b,trackPath, fname);
	unsigned long tPrm  =getFileTime(prmFile);
	unsigned long tTrack=getFileTime(b);
	if(tPrm < tTrack){
		verb("Old profile:  prm_time=%li  track_time=%li\n",tPrm,tTrack); return false;
	}
//======================================================= check parameters
	FILE* f=xopen(prmFile,"rt");
	trackType=getTrackType(fname);
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
			if(strcmp(s2,fname)) {
				failParam("input");  fg=false; break;}}
		else if(strcmp(s1,"bin")==0){
			if(atoi(s2)!=binSize) {failParam("bin");  fg=false; break;}
		}
		else if(strcmp(s1,"bpType")==0){
			if(atoi(s2)!=bpType) {failParam("bpType");  fg=false; break;}
		}
	}
	fclose(f);
	if(chkVersion(bver)){
		failParam("program version"); return false;
	}
	return fg;
}


//=============================================================================
bool bTrack::readPrm(){
	char prmFile[4096], b[4096];
	makeFileName(prmFile,profPath, name, PRM_EXT);
	if(!fileExists(prmFile)) return false;
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
		else if(strcmp(s1,"scale" )==0)  scaleFactor  =atof(s2);
		else if(strcmp(s1,"total" )==0)  total  =atof(s2);
	}
	if(sd==NAN) sd=1;
	fclose(f);
	return true;
}

//=============================================================================
void BuffArray::init(Track *bbt, bool cmpl, bool wrr){		// Prepare
	if(bbt->name==0) return;
	bufBeg=bufEnd=-1; offset=0; wr=false; f=0;
	bt=bbt;

	getMem0(bval,binBufSize+2*wProfSize, "Read bTrack");
	char binFile[4096];
	makeFileName(binFile,profPath, bt->name, BPROF_EXT);
	if(wrr){													// Track for write
		f=xopen(binFile,"w+b"); if(f==0) return;
		for(int i=0; i<binBufSize; i++) bval[i]=NA;			// Fill the file with NA
		if(cmpl) {
			offset=profileLength*sizeof(BINVAL);
			fseek(f,offset,SEEK_SET);
		}
		for(int i=0; i<profileLength; i+=binBufSize){
			int l=profileLength-i; if(l > binBufSize) l=binBufSize;
			fwrite(bval, sizeof(*bval), l, f);
		}
	}
	else{													// Track for read
		f=xopen(binFile,"rb"); if(f==0) return;
		if(cmpl){
			offset=profileLength*sizeof(BINVAL);
			fseek(f,offset,SEEK_SET);
		}
	}
}

void BuffArray::readBuff(int pos){
	bufBeg=(pos/binBufSize)*binBufSize-wProfSize;
	if(bufBeg < 0) bufBeg=0;
	bufEnd=bufBeg+binBufSize+wProfSize*2;
	if(bufEnd > profileLength) bufEnd = profileLength;
	fseek(f,bufBeg*sizeof(BINVAL)+offset,SEEK_SET);
	int l=bufEnd - bufBeg;
	fread(bval, sizeof(BINVAL), l, f);
}

BINVAL BuffArray::get(int pos){
	if(pos < bufBeg || pos >= bufEnd)  {readBuff(pos);}
	BINVAL x=bval[pos-bufBeg];
	return x;
}


void BuffArray::writeBuff(){
	if(wr){		//============= the buffer is changed
		fseek(f,bufBeg*sizeof(*bval)+offset,SEEK_SET);
		int l=bufEnd - bufBeg;
		fwrite(bval, sizeof(*bval), l, f);
	}
}

void BuffArray::set(int pos, BINVAL v){
	if(pos < bufBeg || pos >= bufEnd){
		if(bufBeg >= 0) writeBuff();
		readBuff(pos);
	}
	bval[pos-bufBeg]=v; wr=true;
}

BuffArray::BuffArray(){
	bval=0; f=0; offset=0;	bufBeg=-1; bufEnd=-1; bt=0; wr=false;
}

void BuffArray::close(){
	if(f!=0) {fclose(f);} f=0;
}

BuffArray::~BuffArray(){
	xfree(bval,"~BuffArray");
	close();
}
//==========================================================================
void FloatArray::init(int na){
	if(f==0){                //===================== simple array without file
		for(int i=0; i<profileLength; i++) val[i]=na;	//== fill initial values
		return;
	}
	fseek(f,0,SEEK_SET);
	bufBeg=bufEnd=-1; wr=false;
	for(int i=0; i<binBufSize; i++) val[i]=na;			// Fill the file with NA
	for(int i=0; i<profileLength; i+=binBufSize){
		int l=profileLength-i; if(l > binBufSize) l=binBufSize;
		fwrite(val, sizeof(*val), l, f);
	}
}

FloatArray::~FloatArray(){
	xfree(val,"float profile");
	if(f==0) return;
	fclose(f);
	remove(fname);
	xfree(fname,"~FloatArray file");
}

FloatArray::FloatArray(){
	bufBeg=bufEnd=-1; wr=false; f=0; val=0;
	if(binBufSize>=profileLength){
		getMem0(val,profileLength, "Read bTrack"); return;
	}
	getMem0(val,binBufSize+2*wProfSize, "Read bTrack");
	char binFile[4096];
	unsigned int tt=mtime()+rand();	// generate filename
	sprintf(binFile,"%x.tmp",tt);
	fname=strdup(binFile);
	f=xopen(binFile,"w+b"); if(f==0) return;
}

void FloatArray::readBuf(int pos){
	if(f==0) return;
	//============================================ flush existing buffer before read
	if(bufBeg >= 0)	writeBuf();
	//============================================ define block to read
	bufBeg=(pos/binBufSize)*binBufSize-wProfSize*2;
	if(bufBeg <0) bufBeg=0;
	bufEnd=bufBeg+binBufSize+wProfSize*2;
	if(bufEnd > profileLength) bufEnd = profileLength;
	fseek(f,bufBeg*sizeof(*val),SEEK_SET);
	int l=bufEnd - bufBeg;
	l=fread(val, sizeof(*val), l, f);
	wr=false;
}

float FloatArray::get(int pos){
	//=============================================read the value
	if(pos < bufBeg || pos >= bufEnd)  readBuf(pos);
	float x=val[pos-bufBeg];
	return x;
}
float FloatArray::getLog(int pos){
	//=============================================read the value
	if(pos < bufBeg || pos >= bufEnd)  readBuf(pos);
	float x=val[pos-bufBeg];
	if(x!=NA){
		if(x < 0) x=-log(1-x);
		else      x= log(x+1);
	}
	return x;
}


void FloatArray::writeBuf(){
	if(wr && f){		//=========== file exists & the buffer changed
		fseek(f,bufBeg*sizeof(*val),SEEK_SET);
		int l=bufEnd - bufBeg;
		fwrite(val, sizeof(*val), l, f);
	}
}

void FloatArray::set(int pos, float v){
	if(pos < bufBeg || pos >= bufEnd){	//=== position to write outside the buffer
		readBuf(pos);
	}
	wr=true;
	val[pos-bufBeg]=v;
}
float FloatArray::add(int pos, float v){
	float x=get(pos);
	if(v==NA) return x;
	x=(x==NA) ? v : x+v;
	set(pos,x);
	return x;
}

//=============================================================================
bool bTrack::readBin(){
	//============================================= get file length
	char b[1024];
	if(bytes==0) bytes=new BuffArray();
	makeFileName(b,profPath, name, BPROF_EXT);
	if(!fileExists(b)) return false;
	bytes->init(this,0,0);
	if(hasCompl){
		if(cbytes==0) cbytes=new BuffArray();
		cbytes->init(this,1,0);
	}
	return true;
}

bool Track::openTrack(const char *fname){
	if(!readTrack(fname)) return false;
	getMem0(profWindow,(profWithFlanksLength+10), "bTrack read #2");
	if(doAutoCorr) {
		getMem0(autoCorr,(profWithFlanksLength+10), "bTrack read #2");
		zeroMem(autoCorr, profWithFlanksLength+10);
	}
	return true;
}
bool bTrack::readTrack(const char *fname){
	for(;*fname==DERIV; fname++) {
		deriv++;
	}
	name=strdup(fname);
	//============================================== Read params
	if(!readPrm()) return false;
	if(!readBin()) return false;
	//============================================== Read bytes
	av=sd=nObs=0;
	errStatus=0;
	return true;
}

//===============================================================================
double Track::addStatistics(){
	double avv=0, sdd=0;
//deb(21,"nn=%i",nObs);
	for(int i=0; i<profWithFlanksLength; i++){
		double x=profWindow[i];
		if(x!=NA){avv+=x; sdd+=x*x;  nObs++;}
	}
//deb(22,"nn=%i",nObs);
	av+=avv; sd+=sdd;
	return avv/profWithFlanksLength;
}
void Track::finStatistics(){
	av/=nObs; sd=sd-av*av*nObs; sd/=nObs-1; sd=sqrt(sd);
	if(doAutoCorr){
		double a=autoCorr[0];
		for(int i=0; i<profWithFlanksLength; i++)
			autoCorr[i]/=a;
	}
}

//========================================================================
bool bTrack::isNA(int i, bool cmpl){
	return getBVal(i,cmpl)==NA;
}
bool bTrack::isZero(int i, bool cmpl){
	return getBVal(i,cmpl)<=threshold;
}

void Track::makeIntervals(bool cmpl, IVSet *iv){
	int f=0;
	bool space=false;
	int nZ=0;
	int wp=wProfSize;
	int maxZ=min(maxZero, maxNA);
	for(int i=0; i<wp; i++){
		if(isZero(i,cmpl)) nZ++;
	}
	space = (nZ > maxZ);
	for(int i=wp; i<profileLength; i++){
		if(nZ > maxZ){
			if(!space){
				iv->addIv(f,i-wp);
			}
			space=true;
		}
		else{
			if(space) f=i-wp;
			space=false;
		}
		if( isZero(i   ,cmpl)) nZ++;
		if( isZero(i-wp,cmpl)) nZ--;
	}
	if(!space) {
		iv->addIv(f,profileLength);
	}
}

bool Track::makeIntervals(){
	makeIntervals(0, ivs);
	if(hasCompl) makeIntervals(1, ivsC);
	if(ivs->nIv+ivsC->nIv == 0) {
		writeLog("Track <%s>: no nonZero windows",name);
		fprintf(stderr,"Track <%s>: no nonZero windows",name);
		return false;
	}
	ivs->fin();
	if(ivsC->nIv) ivsC->fin();

	int l0=profileLength; if(hasCompl) l0*=2;
	int l1=ivs->totLength; if(hasCompl) l1+=ivsC->totLength;
	av0=av0*l0/l1;
	sd0=sd0*sqrt(l0/l1);
	return true;
}

//========================================================================
int Track::getRnd(bool cmpl){
	int pos;
	pos=cmpl ? ivsC->randPos() : ivs->randPos();         // get random position in the interval
	return pos;
}


//========================================================================
void Track::clear(){
	ivs->clear();
	if(name) xfree(name,"clear btrack");
}

void bTrack::clear(){
	Track::clear();
	if(bytes ) bytes->close();
	if(cbytes) cbytes->close();
}
//========================================================================
int  Track::countNA(int pos, bool cmpl){
	int c=pos+wProfSize >= profileLength ? wProfSize+pos-profileLength-1 : 0;
	for(int i=0; i< wProfSize && i+pos < profileLength; i++) {
		if(complFg==IGNORE_STRAND){		// ignore strand
			int x=1;
			if(            !isNA(pos+i,0)) x=0;
			if(hasCompl && !isNA(pos+i,1)) x=0;
			c+=x;
		}
		else if(!cmpl || !hasCompl) 		// colinear or profile do not know the orient
			{if(isNA(pos+i,0)) 	c++;}
		else if(cmpl)						// comlement
			{if(isNA(pos+i,1)) 	c++;}
	}
	return c;
}

//========================================================================
int  Track::countZero(int pos, bool cmpl){
	int c=pos+wProfSize >= profileLength ? wProfSize+pos-profileLength-1 : 0;
	for(int i=0; i< wProfSize && i+pos < profileLength; i++) {
		if(            !isNA(pos+i,0) && NAFlag) continue;
		if(hasCompl && !isNA(pos+i,1) && NAFlag) continue;
		if(complFg==IGNORE_STRAND){		// ignore strand
			int x=1;
			if(            !isZero(pos+i,0)) x=0;
			if(hasCompl && !isZero(pos+i,1)) x=0;
				c+=x;
		}
		else if(cmpl || !hasCompl) 		// colinear or profile do not know the orient
			{if(isZero(pos+i,0)) 	c++;}
		else if(cmpl){					// comlement
			{if(isZero(pos+i,1)) 	c++;}
		}
	}
	return c;
}

void Track::init(){
	name=0;
	profWindow=0;
	autoCorr=0;
	avWindow=sdWindow=0;// mean and stdDev in current window
	ivs =new IVSet();
	ivsC=new IVSet();
	av=sd=minP=maxP=nObs=0;
	av0=sd0=total=0;
	deriv=0;
	trackType=0;
	projCoeff=0;
	hasCompl=false;
	scaleFactor=1;
}
void bTrack::initBtr(){
	cbytes=bytes=0;
}

bTrack::bTrack(){
	initBtr();
}
bTrack::bTrack(const char* fname){
	initBtr();
	openTrack(fname);
}

Track::Track(){
	init();
}

Track::~Track(){
	if(profWindow) xfree(profWindow,"~Track 1");
	if(autoCorr)   xfree(autoCorr,"~Track 2");
	if(name)       xfree(name,"~Track 3");
	if(ivs ) del(ivs);
	if(ivsC) del(ivsC);
}
bTrack::~bTrack(){
	if( bytes) del(bytes);
	if(cbytes) del(cbytes);
}
//======================================= decode binary value to real value
double bTrack::getVal(BINVAL b){
	if(b==NA){
			if(NAFlag && trackType==WIG_TRACK){
			double x=rExp();
			double y=x*sd0*noiseLevel;
			return y;
		}
		else return 0;
	}
	if(b <= threshold  && b >= - threshold) return 0;
	double x=b < 0 ? 1-exp(-b/scaleFactor): exp(b/scaleFactor)-1;
	return x;
}

//======================================= decode binary value at certain position
int    bTrack::getBVal(int pos, int cmpl){
	if(pos >= profileLength) return NA;
	if(cmpl){
		if(hasCompl) return cbytes->get(pos);
		else return bytes->get(pos);
	}
	int bv=bytes->get(pos);
	return bv;
}

double bTrack::getValue(int pos, int cmpl){
	BINVAL b=getBVal(pos,cmpl);
	if(b==NA){
			if(NAFlag && trackType==WIG_TRACK){
			double x=rExp();
			double y=x*sd0*noiseLevel;
			return y;
		}
		else return 0;
	}
	if(b <= threshold  && b >= - threshold) return 0;
	double x=exp(b/scaleFactor)-1;
	return x;
}

//====================================== calculate projection to orthogonal subspace
void Track::ortProject(){
	projCoeff=0;
	if(pcorProfile==0) return;
	double xy=0, xx=0;
	for(int i=0; i<profileLength; i++){
		double x,y;
		x=projTrack->getValue(i,false);
		y=getValue(i,false);
		xx+=x*x; xy+=x*y;
		if(hasCompl){
			x=projTrack->getValue(i,true);
			y=getValue(i,true);
			xx+=x*x; xy+=x*y;
		}
	}
	projCoeff=xy/xx;
}

//================================================= Get projected value
double Track::getProjValue(int pos, bool cmpl){
	if(pos <0 || pos >= profileLength) return 0;
	double x=getValue(pos,cmpl);
	if(projCoeff) x-=projCoeff*projTrack->getValue(pos,cmpl);
	return x;
}

//================================================= decode the values to an array
double * Track::getProfile(double *prof, int pos, int l, bool cmpl){ //====== pos - profile position; l -- fragment length without flanks
	int ll=l+LFlankProfSize+RFlankProfSize;
	double *aa=prof+LFlankProfSize;
	double *a=aa;
	avWindow=sdWindow=0;
	//======================================================= fill window
	for(int i=0; i < l; i++, a++){
		double x=0;
		if(complFg==IGNORE_STRAND){				// ignore strand
			x += getProjValue(pos+i,false);		// profile
			if(hasCompl)
				x += getProjValue(pos+i,true);				// + complement profile
		}
		else if(cmpl)                          // colinear or profile do not know the orient
			x += getProjValue(pos+i,hasCompl);
		else
			x += getProjValue(pos+i,false);
		*a=x;

		avWindow+=x; sdWindow+=x*x;
	}

	//======================================= Window Statistics
	avWindow/=l;
	sdWindow=sdWindow-avWindow*avWindow*l;
	//========================================= Constant window => add noise
	a=aa;
	if(sdWindow<=0){
		sdWindow=sd0*1.e-2;
	}
	else{sdWindow/=l; sdWindow=sqrt(sdWindow);}
	//======================================== fill flanks
	int x0=l+LFlankProfSize;
	int x1=x0+LFlankProfSize+RFlankProfSize;

//=========================================	fill flanks
	for(int x=x0; x<x1; x++){
		double xq=rGauss(0,sdWindow)*noiseLevel;
		prof[x%ll]=xq;
	}
	return prof;
}

double * Track::getProfile(int pos, bool cmpl){ //====== pos - profile position; cmpl=true <=> +strand
	return getProfile(profWindow, pos, wProfSize, cmpl);
}
//======================================================
//======================================================
void Track::writeWig(FILE* f, Chromosome *ch){
	int pos1=ch->base;
	int pos2=pos1+ch->length/binSize;
	int fg=0;
	int tr=10;
	for(int i=0,j=pos1; j<pos2; j++,i++){
		double x=getValue(j,0)*4;
		if(x<tr) {fg=0; continue;}
		if(fg==0){
			fprintf(f,"fixedStep chrom=%s start=%i step=%i span=%i\n",ch->chrom,
					i*binSize,binSize,binSize);
		}
		fprintf(f,"%.0f\n",x); fg=1;
	}
}
//======================================================
void Track::writeWig(){
	FILE *f=xopen(name,"wt");
	verb("write track <%s>\n",name);
	fprintf(f,"track type=wiggle_0 description=\"%s\"\n",name);
	for(int i=0; i<n_chrom; i++){
		writeWig(f,chrom_list+i);
	}
	fclose(f);
}

//================================================================================
//================================================================================
//================================================================================
Model::Model(){
	definition=0; form=0; trackType=MODEL_TRACK;
}

//================================================================================
void Model::readModel(const char *fnam){
	name=strdup(fnam);
	char bb[2048];
	char b[2048];
	makeFileName(bb,trackPath, fnam);
	FILE *f=fopen(bb,"r");
	*bb=0;
	for(char *s; (s=fgets(b,sizeof(b),f))!=0;){
		s=trim(s);
		if(*s=='#') continue;
		strcat(bb,b);
	}
	definition=strdup(bb);
	form=frmlInit(bb);
	for(int i=0; i<form->nTracks; i++){
		getTrack(i)->readTrack(getTrackName(i));
	}
}

//================================================================================
bool   Model::readTrack(const char *fname){
	name=strdup(fname);
	readModel(fname);
	return true;
}

//================================================================================
bool   Model::isNA(int pos, bool cmpl){
	for(int i=0; i<form->nTracks; i++){
		bTrack *bt=getTrack(i);
		if(bt->isNA(pos,cmpl)) return true;
	}
	return false;
}
//================================================================================
bool   Model::isZero(int pos, bool cmpl){
	return getValue(pos,0) == 0;
}
//================================================================================
double Model::getValue(int pos, int cmpl){
	return form->calc(pos);
}

//================================================================================
void  Model::clear(){
	for(int i=0; i<form->nTracks; i++){
		bTrack *bt=getTrack(i);
		bt->clear();
	}
}

//================================================================================
Model::~Model(){
	if(definition) xfree(definition,"definition");
	del(form);
}

//================================================================================
//================================================================================
Track * trackFactory(const char* fname){
	if(isModel(fname)){
		Model *mod=new Model();
		mod->openTrack(fname);
		return mod;
	}
	bTrack *bt=new bTrack(fname);
	return bt;
}




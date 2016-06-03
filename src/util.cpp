/*
 * util.cpp
 *
 *  Created on: 17.02.2013
 *      Author: 1
 */

#include "track_util.h"
#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <unistd.h>

#include <sys/file.h>
//#include <dir.h>

const char* version="1.64";

int debugFg=0;
//int debugFg=DEBUG_LOG|DEBUG_PRINT;

const char *debS=0;

Chromosome *chrom_list;       // list of chromosomes
Chromosome *curChrom=chrom_list;
int  stepSize=100;   // frame size fo profile
int  intervFlag0=GENE;   // Flag: for bed tracks uncovered values=0 otherwise uncovered values=NA
bool  NAFlag=0;

long long GenomeLength=0;      // TOTAL LENGTH OF THE GENOME
int n_chrom;
char *trackName;       // current track name
char *chromFile=0;
char *profPath=0;
char *trackPath=0;

char *trackFil=0;		// Track file
char *aliaseFil=0;
char *logFileName=(char*)"./stereogene.log";
char *outFile=0;
char *profile1=0;		// first profile file file name
char *profile2=0;		// second profile file file name
char *resPath=0;
char *statFileName=(char*)"./statistics";
char *paramsFileName=(char*)"./params";
char *mapFil=0;			// map file
char *inputProfiles=0;
char *outTrackFile=0; // Filename for write out track

bool  verbose=0;
bool  silent=0;				// inhibit stdout
bool  syntax=1;				// Strong syntax control

bool  strandFg0=1;
bool  writeDistr=1;
bool  writeBPeak=0;
int  writeDistCorr=TOTAL;		    // write BroadPeak
bool  outSpectr=0;
bool  outChrom=0;
int  outRes=XML|TAB;


int   complFg=0;
int   profileLength;			// size of the profile array
float *profile =0;				// uncompressed profile array
float *profilec=0;				// uncompressed profile array

unsigned char *byteprofile;	 	// compressed profile array
unsigned char *byteprofilec;	// compressed profile array

int   logScale=AUTO_SCALE;

char *pcorProfile=0;    // partial correlation profile file name

int kernelType=KERN_NORM;
double noiseLevel=0.2;
int wSize=100000;        // size of widow (nucleotides)
int wStep=0;             // window step   (nucleotides)
int flankSize=500;
double kernelSigma=1000.;    // kernel width (nucleotides)
double kernelShift=0;      	    // Kernel mean (for Gauss) or Kernel start for exponent
int intervFg0;
double scaleFactor0=0.2;
int outWIG=NONE;
int outThreshold=1;


int wProfStep;          // window step   (profile scale)
int wProfSize;          // size of widow (profile scale)
int LFlankProfSize;        // size of flank (profile scale)
int RFlankProfSize;        // size of flank (profile scale)
int profWithFlanksLength; // size of profWindow array (including random flanks)
double kernelProfSigma;      // kernel width ((profile scale)
double kernelProfShift;
double kernelNS;			// Correction for non-specifisity
bTrack bTrack1, bTrack2, projTrack, mapTrack;
Kernel *kern;
double maxNA0=95;
double maxZero0=95;
double maxNA;
double maxZero;
int nShuffle=100;
int maxShuffle=10000;
int minShuffle=1000;
double pVal=2;
double qVal=0;
MapIv miv;				// interval for mapping
//Map map;
Model model;

int threshold=0;

FILE *logFile=0;
int genFg=0;			// generate data
int lAuto=0;
int corrScale=10;
bool corrOnly=0;
double prod11=0,prod12=0,prod22=0,sprod11=0, sprod12=0,sprod22=0;
Correlation correlation;		// array for correlation picture
Correlation bgcorrelation;		// array for correlation picture
Fourier wCorrelation;

int nProd=0;
int pcaFg=0;
int nPca=100000;
int pcaSegment=100;
double totCorr=0;
unsigned long id;
bool RScriptFg=0;
int bpType=BP_SIGNAL;
int cage=0;
bool clearProfile=false;
int scoreType=AV_SCORE;
AliaseTable alTable;
FileListEntry files[256];
int   nfiles;

double BgAvCorr;
double FgAvCorr;

unsigned int hashx(unsigned int h,unsigned int x);
unsigned int hashx(unsigned int h,char c);
unsigned int hashx(unsigned int h,int x);
unsigned int hashx(unsigned int h,long x);
unsigned int hashx(unsigned int h,const char *s);
unsigned int hashx(unsigned int h,float c);
unsigned int hashx(unsigned int h,double c);

//====================================================================================
ScoredRange::ScoredRange(){
	chr=0; chrom=0;	beg=end=0;	score=0;
}
//====================================================================================
char *AliaseTable::convert(char*oldName){

	char b0[1024],b1[1024];
	strcpy(b0,oldName);

	for(int i=0; i<nAls; i++){
		while(1){
			char* sx=strstr(b0,als[i].oldName);
			if(sx!=0){

				int k=sx-b0;
				char *s0=b0,*s1=b1;
				memcpy(s1,s0,k); s0+=k; s1+=k;
				memcpy(s1,als[i].newName,als[i].lnew);
				s0+=als[i].lold; s1+=als[i].lnew;
				strcpy(s1,s0);

				strcpy(b0,b1);
			}
			else break;
		}
	}
	return strdup(b0);
}

void AliaseTable::readTable(const char* fname){
	FILE *f=gopen(fname,"rt");
	if(f==0) return;
	nAls=0;
	int capacity=100;
	getMem(als,capacity, "AliaseTable::readTable");
	char b[1024],*s;

	for(;(s=fgets(b,sizeof(b),f))!=0;){
		strtok(b,"\r\n");
		s=skipSpace(b);
		if(*s=='#') continue;

		char *s1=strtok(s,"=");
		char *s2=strtok(0," \t");
		s1=strtok(s1," \t");
		if(s1==0 || strlen(s1)==0) continue;
		if(s2!=0) s2=skipSpace(s2);
		else s2=(char*)"";
		als[nAls].oldName=strdup(s1);
		als[nAls].newName=strdup(s2);
		als[nAls].lnew=strlen(s2);
		als[nAls].lold=strlen(s1);
		nAls++;
		if(nAls >= capacity){
			capacity+=100;
			als=(Aliase*)realloc(als,capacity*sizeof(Aliase));
		}
	}
	fclose(f);
}

//====================================================================================
//===============================================  Chromosomes
Chromosome::Chromosome(char *chr, long l, int bb){
    chrom=chr;		//Chromosome name
    length=l;		//Chromosome length
    base=bb;			//Start position in binary profile
    av1=av2=corr=lCorr=count=0;
    densCount=0;
    distDens=0;
}

void Chromosome::clear(){
    av1=av2=corr=lCorr=count=0; densCount=0;
    if(profWithFlanksLength) {
    	getMem0(distDens,profWithFlanksLength,  "Chromosome::clear #1");
    	zeroMem(distDens,profWithFlanksLength);
    }
}

void clearChromosomes(){
	for(int i=0; i<n_chrom; i++){
		chrom_list[i].clear();
	}
}


int CHROM_BUFF=30;

int readChromSizes(char *fname){

	if(fname==0) errorExit("Chromosome file undefined");
	FILE *f=xopen(fname,"rt");
	if(f==0)return 0;
	char buff[2048];
	n_chrom=0;
	int max_chrom=CHROM_BUFF;
	chrom_list=(Chromosome *)malloc(max_chrom*sizeof(Chromosome));
	char *s1,*s2;

	for(char *s; (s=fgets(buff,2048,f))!=0;){
		if (strcmp(s, "") == 0) continue;
		if(n_chrom>=max_chrom) {
			max_chrom+=CHROM_BUFF;
			chrom_list=(Chromosome *)realloc(chrom_list,max_chrom*sizeof(Chromosome));
		}
	    if((s1=strtok(s," \t\r\n"))==NULL) continue;
	    if((s2=strtok(0," \t\r\n"))==NULL) continue;

	    chrom_list[n_chrom]=Chromosome(strdup(s1), atol(s2), profileLength);

		int filLen=(chrom_list[n_chrom].length+stepSize-1)/stepSize;
		profileLength+=filLen;

		GenomeLength+=chrom_list[n_chrom].length;
	    n_chrom++;

	}
	curChrom=chrom_list;
	return 1;
};

Chromosome* findChrom(char *ch){
	if(strcmp(curChrom->chrom, ch) ==0) return curChrom;
	for(int i=0; i<n_chrom; i++){
		curChrom=chrom_list+i;
		if(strcmp(curChrom->chrom, ch) ==0) return curChrom;
	}
	fprintf(stderr,"Chromosome %s not found\n",ch);
	return 0;
}
//========================================================================================
long pos2filePos(char*chrom,long pos){
	Chromosome *ch=findChrom(chrom);
	if(ch==0) return 0;
	curChrom=ch;
	long p=pos/stepSize+curChrom->base;
	return p;
}

//========================================================================================
Chromosome *getChromByPos(int pos){
	Chromosome *ch0=chrom_list;
	for(int i=0; i<n_chrom; i++){
		if(pos < chrom_list[i].base) return ch0;
		ch0=chrom_list+i;
	}
	return ch0;
}
void filePos2Pos(int pos, ScoredRange *gr, int length){

	Chromosome *ch0=getChromByPos(pos);


	if(ch0==0) return;
	pos-=ch0->base;
	long p1=stepSize*long(pos);
	gr->chr=ch0;
	gr->chrom=ch0->chrom;
	gr->end=(gr->beg=p1)+length;
	return;
}
//========================================================================================
Chromosome *checkRange(ScoredRange *gr){
	Chromosome* chr=findChrom(gr->chrom);
	if(chr==0) return 0;
	if((gr->beg < 0) || (gr->end >= chr->length)) return 0;
	return chr;
}
//========================================================================================
char *skipSpace(char *s)  {while(*s!=0 &&  isspace(*s)) s++; return s;}
char *skipNoSpace(char *s){while(*s!=0 && !isspace(*s)) s++; return s;}
int EmptyString(const char*buff){
	for(const char*s=buff; *s!=0; s++)
		if(!isspace(*s)) return 0;
	return 1;
}

//==============================================================================
bool isUInt(const char *s){
	for(;*s;s++) {
		if(!isdigit(*s)) return false;
	}
	return true;
}
bool isInt(const char *s){
	char *ss=skipSpace((char*)s);
	if(*ss=='-' || *ss=='+') ss++;
	return isUInt(ss);
}
bool isDouble(const char *s){
	char b[256]; strcpy(b,s);
	char *s0=strtok(b,"eE");
	char *s1=strtok(0,"");
	if(s1!=0 && strlen(s1)!=0 && !isInt(s1)) return false;
	s0=strtok(s0,"."); s1=strtok(0,".");
	if(!isInt(s0)) return false;
	if(s1!=0 && strlen(s1)!=0 && !isUInt(s1)) return false;
	return true;
}

// extract attribute value by attr name
char * getAttr(char *s0, const char *name, char *buf){
	char *s=s0;
	while(*s!=0){
		char *ss=buf;
		s=strstr(s,name);
		if(s==0) return 0;
		s=skipSpace(s+strlen(name));
		if(*s==0 || *s!='=') continue;
	    s=skipSpace(s+1);
		if(*s!='\"'){while(*s!=0 && !isspace(*s)) *ss++=*s++;}
		else{   s++; while(*s!=0 && *s != '\"'  ) *ss++=*s++;}
		*ss=0; return buf;
	}
	return 0;
}

// convert string to upper case
char *strtoupper(char*s){
	for(char *ss=s;*ss;ss++) *ss=toupper(*ss);
	return s;
}
// Errors
//================
const char *errStatus=0;

FILE *openLog(){
	if(logFileName) {
		FILE *f=gopen(logFileName,"at");
		if(f!=0) return f;
		else{
			fprintf(stderr, "Error in opening log file%s\n", logFileName);
		}
	}
	return 0;
}

void errorExit(const char *format, va_list args){
    fflush(stdout);
	if (format != NULL) {
		char b[1024];
		vsprintf(b, format, args);
	    fprintf(stderr, "%s", b);
	    if(errStatus) fprintf(stderr, "%s\n", errStatus);
	    else fprintf(stderr, "\n");
	    FILE *f=openLog();
		if(f) {
			flockFile(f);
			fprintf(f,"#%8lx-> %s",id,b);
			if(errStatus) fprintf(f, "%s\n", errStatus);
			else fprintf(f, "\n");
			funlockFile(f);
			fclose(f);
		}
	}
	exit(-1);
}
void errorExit(const char *format, ...){
	va_list args;
	va_start(args, format);
	errorExit(format, args);
	va_end(args);
}

void writeLog(const char *format, va_list args){
	FILE *f=openLog();
	if(f) {
		flockFile(f);
		fprintf(f,"#%8lx-> ",id);
		vfprintf(f,format,args);
		funlockFile(f);
		fclose(f);
	}
}

void writeLog(const char *format, ...){
	va_list args;
	va_start(args, format);
	writeLog(format, args);
	va_end(args);
}

void verb_(const char *format, va_list args){
	if(verbose){
		vprintf(format,args);
	}
}
void verb(const char *format, ...){
	va_list args;
	va_start(args, format);
	verb_(format, args);
	va_end(args);
	fflush(stdout);
}
void xverb_(const char *format, va_list args){
	if(!silent){
		vprintf(format,args);
	}
}
void xverb(const char *format, ...){
	va_list args;
	va_start(args, format);
	xverb_(format, args);
	va_end(args);
	fflush(stdout);
}

//============================================ Debugging
//======================================================
FILE *debLogFile=0;
void clearDeb(){
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","wt");
		fclose(debLogFile); debLogFile=0;
	}
}
//===========================================
void _deb_(const char *format, va_list args){
	char b[2048];
	vsprintf(b, format, args);
	if((debugFg&DEBUG_PRINT)!=0){
		printf("%s\n",b);
	}
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","a+t");
		fprintf(debLogFile, "%s\n",b);
		fclose(debLogFile); debLogFile=0;
	}
}
//===========================================

void deb_(int num){
	if((debugFg&DEBUG_PRINT)!=0){
		if(debS) printf("%s",debS);
		printf(" #%i ",num);
	}
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","a+t");
		if(debS) fprintf(debLogFile, "%s",debS);
		fprintf(debLogFile, " #%i ",num);
		fclose(debLogFile); debLogFile=0;
	}
}

void deb(int num){
	if((debugFg&DEBUG_PRINT)!=0){
		if(debS) printf("%s",debS);
		printf(" #%i\n",num);
	}
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","a+t");
		if(debS) fprintf(debLogFile,"%s", debS);
		fprintf(debLogFile, " #%i\n",num);
		fclose(debLogFile); debLogFile=0;
	}
}


void deb(const char *format, ...){
	va_list args;
	va_start(args, format);
	_deb_(format, args);
	va_end(args);
}

void deb(int num, const char *format, ...){
	if(debugFg==0) return;
	deb_(num);
	va_list args;
	va_start(args, format);
	_deb_(format, args);
	va_end(args);
}
//====================================================================================
char *Timer::getTime(){
	long dt=getTimer();
	int ms=(int)(dt%1000);
	int s=(int)(dt/1000);
	sprintf(bb,"%i.%is",s,ms);
	return bb;
}

Timer::Timer(){
	reset();
}
void Timer::reset(){
	start=mtime();
}

long Timer::getTimer(){
	long curTime=mtime();
	return curTime-start;
}

long Timer::mtime()
{
  struct timeval t;

  gettimeofday(&t, NULL);
  long mt = (long)t.tv_sec * 1000 + t.tv_usec / 1000;
  return mt;
}


//======================================================

// open file with control
//================
char* parseTilda(char *b, const char*fname){
	if(*fname=='~'){
		char *z=getenv("HOME");
		if(z==0) z=getenv("HOMEPATH");
		if(z!=0) {
			strcpy(b,correctFname(z));
			if(b[strlen(b)-1] != '/') strcat(b,"/");
			fname+=2;}
	}
	else *b=0;
	return strcat(b,fname);
}

FILE *xopen(const char* fname, const char *t){
	if(fname==0) errorExit("can\'t open file <null>");
	FILE* f=gopen(fname,t);
	if(f==0){
		char b[2048];
		errorExit("can\'t open file <%s> (<%s>)",fname, parseTilda(b,fname));
	}
	return f;
}

FILE *gopen(const char*fname, const char* type){		// open file with parsing ~
	char b[2048];
	return fopen(parseTilda(b,fname),type);
}

//===================

void *xmalloc(size_t n, const char *err){
	void *a=malloc(n);
	if(a==0){
		if(err==0)
			errorExit("can't allocate memory: %li",n);
		else
			errorExit("can't allocate memory in %s: %li",err,n);
	}
	return a;
}
//===================
// standard random gaussian variable
double rGauss(){
	double phi=(double)rand()/RAND_MAX * 2 * M_PI, r=0;
	while(r==0) r=(double)rand()/RAND_MAX;
	r=sqrt(-2.*log(r));
	return r*sin(phi);
}

// random gaussian variable with given mean anf std deviation
double rGauss(double e, double sigma){
	return rGauss()*sigma+e;
}

// random integer in given interval
unsigned long randInt(unsigned long n){
	unsigned long rn=(unsigned long)((double)(rand())/(double(RAND_MAX))*n);
	return rn;
}

// remove fucked backslash
char *correctFname(const char* s){
	char *b=strdup(s),*ss;
	for(ss=b;*ss;ss++) if(*ss=='\\') *ss='/';
	return b;
}

Histogram::Histogram(int n){
	minVal=-1; maxVal=1; e=0; sigma=0; beta=1; iq=iqq=im=0;
	nBin=n; count=0;
	dd=0; db=0; Fp=0; Fm=0;
	errStatus="init Histogram";
	getMem(dd, nBin, "Histogram #1");
	getMem(db, nBin, "Histogram #1");
	getMem(Fp, (nBin+1), "Histogram #1");
	getMem(Fm, (nBin+1), "Histogram #1");
	zeroMem(dd,nBin);
	bin=(maxVal-minVal)/nBin;
	ready=false;
	errStatus=0;
}
void Histogram::add(double x){
	if(ready) return;
	int i=int((x-minVal)/bin);
	if(i<0) i=0; if(i>=nBin) i=nBin-1;
	dd[i]++; e+=x; sigma+=x*x;
	count++;
}

void Histogram::normBeta(){		// Calculate cummulative Beta distribution
//	printf("normBeta %i\n",ready);
	if(ready) return; ready=true;
	norm();
	fitBeta();
	double eb=0;
	for(int i=0; i<nBin; i++){
		double x=minVal+bin*i;
		eb+=(db[i]=xBetaD(beta,beta,x));
	}
	for(int i=0; i<nBin; i++){db[i]/=eb*bin;}
	Fp[0]=Fm[nBin]=0; Fm[0]=1;
	for(int i=1, j=nBin-1; i<nBin; i++,j--){
		Fp[i]=Fp[i-1]+db[i]*bin;
		Fm[j]=Fm[j+1]+db[j]*bin;
	}

}
void Histogram::normF(){
	if(ready) return; ready=true;
	norm();
	Fp[0]=Fm[nBin]=0; Fm[0]=1;
	for(int i=1, j=nBin-1; i<nBin; i++,j--){
		Fp[i]=Fp[i-1]+dd[i]*bin;
		Fm[j]=Fm[j+1]+dd[j]*bin;
	}
}

void Histogram::norm(){
	e/=count; sigma=(sigma-count*e*e)/(count-1);
	beta=(1/sigma-1)/2;
	sigma=sqrt(sigma);
	for(int i=0; i<nBin; i++){
		dd[i]=dd[i]/count/bin;
	}
}

double Histogram::interpol(double x, double *fun){
	int i=int((x-minVal)/bin);
	if(i<0) {i=0;}
	if(i>=nBin) {i=nBin-1;}
	double dx0=x-(minVal+bin*i), dx1=bin-dx0;
	double f0=log(fun[i]), f1=log(fun[i+1]);
	return exp((f0*dx1+f1*dx0)/bin);
}
//=======================================================
double Histogram::error(double b){
	double z=0;
	for(int i=nBin-1; i>iqq; i--){
		double x=bin*i + minVal;
		z+=pow(1-x*x,b)*bin;
	}
	z=0.25/z;
	double e=0,s=0;
	for(int i=nBin-1; i>iqq; i--){
		double x=bin*i + minVal;
		x=z*pow(1-x*x,b)-dd[i];
		e+=x;
		s+=x*x;
	}
	s*=bin; s=sqrt(s);
	return s;
}


void Histogram::fitBeta(){
	double s=0;
	for(int i=0; i<nBin; i++){
		s+=dd[i]*bin;
		if(s<0.25) iq=i;
		if(s<0.5 ) im=i;
		if(s<0.75) iqq=i;
	}
//	beta=1;
	double d=1;
	double d1=error(beta-d);
	double d2=error(beta);
	double d3=error(beta+d);
	int k=0;
	while(d > 0.01 && k<1000){
		if(d1<=d2) {
			beta-=d;
			d3=d2; d2=d1; d1=error(beta-d);
		}
		else if(d3<=d2){
			beta+=d;
			d1=d2; d2=d3; d3=error(beta+d);
		}
		else{
			double nn=(d3-d1), dd=d3+d1-2*d2;
			if(dd==0) dd=1;
			dd=d*nn/dd; if(dd>2) dd=2; if(dd<-2) dd=-2;
			beta-=dd; d=abs(dd);
			d1=error(beta-d); d2=error(beta); d3=error(beta+d);
			if(d2<=d1 && d2<=d3) d/=3;
			else d*=2;
			if(d>2) d=randInt(4);
		}
//deb("d=%f  k=%i",d,k);
		k++;
	}
}

double Histogram::pValp(double x){return interpol(x,Fp);}
double Histogram::pValm(double x){return interpol(x,Fm);}

double xBetaD(double a, double b, double x){
	return pow(1-x,a)*pow(1+x,b);
}

void Histogram::print(FILE *f){
	for(int i=0; i<nBin; i++){
		double x=minVal+bin*i;
		fprintf(f,"%f\t%f\t%.2e\t%.2e\t%.2e\n",x,dd[i], db[i],Fp[i],Fm[i]);
	}
}

double norm(double *x, int l){
	double d=0,e=0,dd,ee;
	for(int i=0; i<l; i++){d+=x[i]*x[i]; e+=x[i];}
	ee=e/l; d=d*l-e*e; dd=d/((l-1)*l);
	if(dd<0) dd=0; dd=sqrt(dd);
	if(dd <= ee*ee*1.e-5) {
		return 0;}

	for(int i=0; i<l; i++) x[i]=(x[i]-ee)/dd;
	return dd;
}

FileName::FileName(char *fn){
	char bb[1024]; strcpy(bb,fn);
	char *x=strrchr(bb,'/');
	if(x==0) {x=bb; path=".";}
	else {*x=0; x++; path=strdup(bb);}
	char *s=strrchr(x,'.');
	if(s==0){ext="";}
	else	{*s=0; s++; ext=strdup(s);}
	name=strdup(x);
}

char *FileName::fname(){
	char bb[1024]; sprintf(bb,"%s/%s.%s",path,name,ext);
	return strdup(bb);
}
//================ create filename without extension
char *FileName::coreFname(){
	char bb[1024]; sprintf(bb,"%s/%s",path,name);
	return strdup(bb);
}

//================= create filename using path and name
char* makeFileName(char *b, const char *path, const char*fname){
	if(path==0) return strcpy(b,fname);
	if(*fname=='/' || *fname=='~') return strcpy(b,fname);
	char *s;
	if((s=strrchr((char*)fname,'/'))!=0) fname=s+1;
	sprintf(b,"%s%s",path,fname);
	return b;
}
char *makeFileName(char *b, const char *path, const char*fname, const char*ext){
	makeFileName(b,path,fname);
	char *ss=strrchr(b,'/'); if(ss==0) ss=b;
	char *s=strrchr(ss,'.'); if(s) *s=0;
	return strcat(strcat(b,"."),ext);
}

int _makeDir(const char * path){
    struct stat sb;
    if (stat(path, &sb) == 0 && S_ISDIR(sb.st_mode)) return 0;
#if defined(_WIN32)
	return _mkdir(path);
#else
	mode_t mode=S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
	return mkdir(path, mode); // notice that 777 is different than 0777
#endif
}
void makeDir(const char *path){
	char b[2048];
	parseTilda(b,path);
	char *s=b+strlen(b)-1;
	if(*s=='/') *s=0;
	for(char *s=b; (s=strchr(s+1,'/'))!=0;){
		*s=0;
		if(_makeDir(b))
			errorExit("Can not create directory %s\n",b);
		*s='/';
	}
	_makeDir(b);
}
//================== make path - add '/' if necessary
char* makePath(char* pt){
	if(pt==0) return pt;
	char b[2048];
	char *s=pt+strlen(pt)-1;
	if(*s=='/') *s=0;
	return strdup(strcat(strcpy(b,pt),"/"));
}

//=================== change file extension
char * cfgName(char* p, char* ext){
	FileName fn(p);	fn.ext=ext;
	return fn.fname();
}
void flockFile(FILE *f){
#if defined(_WIN32)
	return;
#else
	flockfile(f);
#endif
}
void funlockFile(FILE *f){
#if defined(_WIN32)
	return;
#else
	funlockfile(f);
#endif
}


//=================== Check if given file exists
bool fileExists(const char *fname){
	bool fg=false;						// The file do not exist. The header should be writen.
	FILE *f=gopen(fname,"rt");	// check if statistics file exists
	if(f!=0) {fg=true; fclose(f);}
	return fg;
}

bool fileExists(const char* path, const char *fname){				// check if the file exists
	char b[4096];
	makeFileName(b,path,fname);
	return fileExists(b);
}
bool fileExists(const char* path, const char *fname, const char *ext){ // check if the file exists
	char b[4096];
	makeFileName(b,path,fname,ext);
	return fileExists(b);
}
unsigned long getFileTime(const char *fname){
	struct stat mystat;
	stat(fname,&mystat);
	return mystat.st_mtime;
}

const char *getIvFlag(){
				switch(intervFlag0){
				case GENE:     return "gene"; break;
				case EXON:     return "exon"; break;
				case IVS:      return "ivs" ; break;
				case GENE_BEG: return "gene_beg"; break;
				case EXON_BEG: return "exon_beg"; break;
				case IVS_BEG:  return "ivs_beg" ; break;
				case GENE_END: return "gene_end"; break;
				case EXON_END: return "exon_end"; break;
				case IVS_END:  return "ivs_end" ; break;
				default: return ""; break;
				}
				return "";
}
//======================================================================
int nearPow2(int n, int &i){
	int nn=1;
	for(i=0; i<30; i++){
		if(nn >= n) return nn;
		nn*=2;
	}
	return nn;
}
int nearPow2(int n){
	int z;
	return nearPow2(n,z);
}

int nearFactor(int n){
	int qMin;
	int pow2=nearPow2(n);
	qMin=pow2;
	int k3=1,k35=1;
	for(k3=1; k3<pow2; k3*=3){
		for(k35=k3; k35 < pow2; k35*=5){
				int p2=nearPow2(k35);
				int qq=pow2*k35/p2;
				if(qq>=n && qq<qMin) qMin=qq;
		}
	}
return qMin;
}

//=================== extract file name
char *getFname(char *s){
	char *ss=strrchr(s,'/');
	if(ss==0) ss=s; else ss++;
	return ss;
}

const char *getExt(const char *fname){
	const char *s=strrchr(fname,'/');
	if(s==0) s=fname;
	s=strrchr(s,'.');
	if(s==0) return 0;
	return s+1;
}
char *getFnameWithoutExt(char *buf, char *fname){
	char *s;
	s=strrchr(fname,'/'); if(s==0) s=fname; else s++;
	strcpy(buf,s);

	s=strrchr(buf,'.'); if(s) *s=0;
	return buf;
}

//========================================================================================
// search appropriate cfg file
//void readCfg(int argc, const char *argv[]) {
//	argv[0]=correctFname(argv[0]);
//	char *cfg=cfgName((char*)argv[0], (char*)"cfg");
//	readCfg(cfg);					// deafult cfg
//	char* cfg1=strrchr(cfg,'/');	// cfg in current directory
//	if(cfg1 !=0) readCfg(cfg1+1);
//	for(int i=0; i<argc; i++){
//		if(strncmp(argv[i],"cfg=",4)==0) {
//			verb("read cfg <%s>\n",cfg);
//			readCfg((char*)(argv[i]+4));
//		}
//	}
//}

int getFlag(char*s){
	int fg=0;
	if(		keyCmp(s,"1")==0 || keyCmp(s,"YES")==0 || keyCmp(s,"ON" )==0) {fg=1;}
	else if(keyCmp(s,"0")==0 || keyCmp(s,"NO")==0  || keyCmp(s,"OFF")==0) {fg=0;}
	else fg=-1;
	return fg;
}

int keyCmp(const char *str, const char *key){
	for(;;str++, key++){
		if(*str==0){
			if(*key==0) return 0;
			else return 1;
		}
		if(*key==0) return -1;
		if(toupper(*str) != toupper(*key)) return 1;
	}
	return strncmp(str,key,strlen(key));
}

const char*getKernelType(){
	const char *type;
	if(kernelType==KERN_NORM	 ) type="N";
	else if(kernelType==KERN_LEFT_EXP ) type="L";
	else if(kernelType==KERN_RIGHT_EXP) type="R";
	else return "X";
	char b[80];
	if(kernelShift >0){
		sprintf(b,"%s_%.1fR",type,kernelShift/1000);
		return strdup(b);
	}
	else if(kernelShift <0){
		sprintf(b,"%s_%.1fL",type,-kernelShift/1000);
		return strdup(b);
	}
	else return type;

}

int   fileId=0;

void addFile(const char* fname, int id){
	files[nfiles].fname=strdup(fname);
	files[nfiles].id=fileId;
	nfiles++;
}

void addFile(const char* fname){
	if(nfiles > 256) errorExit("too many input files\n");
	char b[4096], *s;
	strcpy(b,fname); s=strrchr(b,'.'); if(s) s++;
	if(s && (keyCmp(s,"lst")==0 || keyCmp(s,"list")==0)){
		FILE *f=0;
		if(fileExists(fname)) f=xopen(fname,"rt");
		else{
			makeFileName(b,trackPath,fname);
			f=xopen(b,"rt");
		}
		for(;(s=fgets(b,sizeof(b),f))!=0;){
			strtok(b,"\r\n#");
			s=skipSpace(b);
			if(strlen(s)==0 || *s=='#') continue;
			addFile(s, fileId);
		}
		fclose(f); fileId++;
		return;
	}
	else{
		addFile(fname, fileId++);
	}
}


char *trim(char *s){
	if(s==0) return 0;
	s=skipSpace(s);
	for(int i=strlen(s)-1; i>=0; i--){
		if(isspace(s[i])) s[i]=0;
		else break;
	}
	if(*s=='\"') s++;
	char *ss=strchr(s,'\"'); if(ss) *ss=0;
	return s;
}


unsigned int hashx(unsigned int h,char c){
	return h+(c-32)+1234567;
}
unsigned int hashx(unsigned int h,const char *s){
	if(s==0) return h;
	for(;*s;s++) h=hashx(h,*s);
	return h;
}
unsigned int hashx(unsigned int h,unsigned int x){
	return h*3+x;
}
unsigned int hashx(unsigned int h,int x){
	return h*3+x;
}
unsigned int hashx(unsigned int h,long x){
	return h*3+x;
}
unsigned int hashx(unsigned int h,float c){
	unsigned int f; memcpy(&f,&c,sizeof(float));
	return hashx(h,f);
}
unsigned int hashx(unsigned int h,double c){
	float f=(float) c;
	return hashx(h,f);
}

void makeId(){
	id=hashx(id,outFile);
	id=hashx(id,chromFile);
	id=hashx(id,stepSize);
	id=hashx(id,intervFlag0);
	id=hashx(id,qVal);
	id=hashx(id,pVal);
	id=hashx(id,nShuffle);
	id=hashx(id,maxZero);
	id=hashx(id,maxNA0);
	id=hashx(id,kernelSigma);
	id=hashx(id,profile1);
	id=hashx(id,profile2);
	id=hashx(id,wSize);
	id=hashx(id,wStep);
	id=hashx(id,flankSize);
	id=hashx(id,kernelType);
}


//================= Create directories
void makeDirs(){
	if(profPath!=0) makeDir(profPath);
	else profPath=strdup("./");
	if(resPath!=0) makeDir(resPath);
	else resPath=strdup("./");
	if(trackPath!=0) makeDir(trackPath);
	else trackPath=strdup("./");
}


//int main0(int argc, const char *argv[]) {
////	clearDeb();
////	debugFg=DEBUG_LOG|DEBUG_PRINT;
//
//	char *chrom=getenv("SG_CHROM");
//	if(chrom!=0) chromFile=strdup(chrom);
//
//	unsigned long t=time(0);	id=t&0xffffff;	// define run id
//	readArgs( argc, argv);
//	if(wStep==0)   wStep=wSize;
//	if(RScriptFg) {writeDistCorr|=TOTAL; writeDistr=1;}
//	if(complFg==0){
//		if(strandFg0) complFg=COLLINEAR;
//		else		  complFg=IGNORE_STRAND;
//	}
//
//
//
//	FileName fn((char*) argv[0]);
//	const char *progName=fn.name;
//
//	makeDirs();
//
//	verb("===== %s version %s =====\n",progName,version);
//	if(nfiles==0){
//		printf("\n");
//		printf("Program %s compares two track and calculates Kernel Correlation\n",progName);
//		printf("Usage:\n");
//		printf("$ ./%s [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n",progName);
//		printf("\n");
//		exit(0);
//	}
//
//	if(aliaseFil!=0)  alTable.readTable(aliaseFil);		// read aliases
//	readChromSizes(chromFile);							// read chromosomes
//
//	if(cage) {CageMin(files[0].fname,files[1].fname); exit(0);}
//
//	if(pcaFg) pcaMain(profile1);
//
//	Correlator();
//
//	if(logFile) {fclose(logFile); logFile=0;}	// debug log file
//	return 0;
//}

/*
 * util.cpp
 *
 *  Created on: 17.02.2013
 *      Author: Mironov
 */
#include "track_util.h"
#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <unistd.h>


#include <sys/file.h>


char *logFileName=(char*)"./stereogene.log";
unsigned long id;
//char *outPath=0;
char *repFile=0;		// output filename
bool  verbose=0;
bool  silent=0;				// inhibit stdout
//int debugFg=0;
int debugFg=0;
const char *debS=0;


//======================================================================================
//============================    String  Parsing ======================================
//======================================================================================
char *skipSpace(char *s)  {while(*s!=0 &&  isspace(*s)) s++; return s;}
const char *skipNoSpace(const char *s){while(*s!=0 && !isspace(*s)) s++; return s;}
int isEmpty(const char*buff){
	for(const char*s=buff; *s!=0; s++)
		if(!isspace(*s)) return 0;
	return 1;
}


//==============================================================================
//==============================================================================
char *lastChar(char *s){
	if(s==0 ||*s==0) return 0;
	return s+strlen(s)-1;
}
//==============================================================================
bool isUInt(const char *s){
	for(;*s;s++) {
		if(!isdigit(*s)) break;
	}
	if(*s==0) return true;
	if(s[1]==0){
		if(toupper(*s)=='K') return true;
		if(toupper(*s)=='M') return true;
	}
	return false;
}
//==============================================================================
bool isInt(const char *s){
	char *ss=skipSpace((char*)s);
	if(*ss=='-' || *ss=='+') ss++;
	return isUInt(ss);
}
//==============================================================================
bool isDouble(const char *s){
	char b[256]; strcpy(b,s);
	char *last=lastChar(b);
	char c=toupper(*last);
	if(c=='K' || c=='M') *last=0;

	char *s0=strtok(b,"eE");
	char *s1=strtok(0,"");
	if(s1!=0 && strlen(s1)!=0 && !isInt(s1)) return false;
	s0=strtok(s0,"."); s1=strtok(0,".");
	if(!isInt(s0)) return false;
	if(s1!=0 && strlen(s1)!=0 && !isUInt(s1)) return false;
	return true;
}

int readInt(const char* val){

	return int(readDouble(val));
}

double readDouble(const char *s){
	char b[256]; strcpy(b,s);
	double k=1;
	char *last=lastChar(b);
	char c=toupper(*last);
	if(c=='K') {k=1000; *last=0;}
	if(c=='M') {k=1000000; *last=0;}

	double v=atof(b)*k;
	return v;
}


//=================================== extract attribute value by attr name
char * getAttr(char *s0, const char *name, char *buf){
	char *s=s0;
	int ll=strlen(name);
	char *ss=buf;
	while(*s!=0){
		s=strchr(s,*name);
		if(s==0) return 0;
		if(strncmp(s,name,ll)!=0) { s++; continue;}
		s=skipSpace(s+strlen(name));
		if(*s==0 || *s!='=') continue;
	    s=skipSpace(s+1);
		if(*s!='\"'){while(*s!=0 && !isspace(*s)) *ss++=*s++;}
		else{   s++; while(*s!=0 && *s != '\"'  ) *ss++=*s++;}
		*ss=0; return buf;
	}
	return 0;
}


//=================================== convert string to upper case
char *strtoupper(char*s){
	for(char *ss=s;*ss;ss++) *ss=toupper(*ss);
	return s;
}


//===================== check if given string contains given key: 0 -- contains; 1 -- does not
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
//========================= trim given string
char *trim(char *s){
	if(s==0) return 0;
	s=skipSpace(s);
	for(int i=strlen(s)-1; i>=0; i--){
		if(isspace(s[i])) s[i]=0;
		else break;
	}
	if(*s=='\"') s++;
	char *ss=strrchr(s,'\"'); if(ss) *ss=0;
	return s;
}

char* skipInt(char *s){
	if(*s=='-' || *s=='+') s++;
	while(isdigit(*s)) s++;
	return s;
}


bool isfloat(char *s){
	s=skipInt(s);
	if(*s==0) return true;
	if(*s=='.') s=skipInt(s+1);
	if(*s==0) return true;
	if(*s != 'e' && *s !='E') return false;
	s=skipInt(s+1);
	if(*s==0) return true;
	return false;
}

//==================================================================
Alias::Alias(const char *o, const char *r){
	old=strdup(o);
	rpl=strdup(r);
	l_rpl=strlen(rpl);
	l_old=strlen(old);
}
int Alias::replace(char *txt){
	char *s;
	for(int n=0;;n++){
		s=strstr(txt,old);
		if(s==0) return n;
		char b[4096];
		strcpy(b,s+l_old);
		strcpy(s,rpl);
		strcat(s,b);
		txt=s+l_rpl;
	}
	return 0;
}

//============================================================
//=======================    Logging    ======================
//============================================================
const char *errStatus=0;
void clearLog(){
	if(logFileName) fclose(xopen(logFileName,"wt"));
}


FILE *openLog(){
	if(logFileName) {
		FILE *f=xopen(logFileName,"at");
		if(f!=0) return f;
		else{
			fprintf(stderr, "Error in opening log file%s Error code=%i\n", logFileName, errno);
		}
	}
	return 0;
}


void writeLog(const char *format, va_list args){
	FILE *f=openLog();
	if(f) {
		flockFile(f);
		if(!debugFg) fprintf(f,"#%08lx-> ",id);
		fprintf(f,"#%i-> ",getpid());
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
void writeLogErr(const char *format, ...){
	char b[TBS];
	va_list args;
	va_start(args, format);
	vsprintf(b, format,args);
	va_end(args);
	writeLog(b);
	fprintf(stderr,"%s",b);
}


//============================================================
//=======================    Verbose    ======================
//============================================================
void verb_(const char *format, va_list args){	//======== write if verbose =1
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
void xverb_(const char *format, va_list args){   //======== write if silent =0
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


//============================================================
//=======================    Errors    =======================
//============================================================
void errorExit(const char *format, va_list args){
    fflush(stdout);
	if (format != NULL) {
		char b[1024];
		vsprintf(b, format, args);
	    fprintf(stderr, "%s", b);
	    if(errStatus) fprintf(stderr, " %s\n", errStatus);
	    else fprintf(stderr, "\n");
	    FILE *f=openLog();
		if(f) {
			flockFile(f);
			fprintf(f,"#%08lx-> %s",id,b);
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




//============================================================
//=======================    Debug     =======================
//============================================================
Timer debTimer;
FILE *debLogFile=0;
void clearDeb(){
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","wt");
		fclose(debLogFile); debLogFile=0;
	}
}
//===========================================
void _deb_(bool t, const char *format, va_list args){
	char b[TBS];
	vsprintf(b, format, args);
	if((debugFg&DEBUG_PRINT)!=0){
		printf("%s",b);
		if(t) printf(" %s",debTimer.getTime());
		printf("\n");
	}
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","a+t");
		fprintf(debLogFile, "%s",b);
		if(t) fprintf(debLogFile, " %s",debTimer.getTime());
		fprintf(debLogFile, "\n");
		fclose(debLogFile); debLogFile=0;
	}
	if(t) debTimer.reset();
}
//===========================================


void deb(int num, bool t, char e){
	if((debugFg&DEBUG_PRINT)!=0){
		if(debS) printf("%s",debS);
		printf(" #%i",num);
		if(t) printf(" %s",debTimer.getTime());
		printf("%c",e);
	}
	if((debugFg&DEBUG_LOG)!=0){
		if(debLogFile==0) debLogFile=fopen("deb_log","a+t");
		if(debS) fprintf(debLogFile, "%s",debS);
		fprintf(debLogFile, " #%i",num);
		if(t) fprintf(debLogFile, " %s",debTimer.getTime());
		fprintf(debLogFile, "%c",e);
		fclose(debLogFile); debLogFile=0;
	}
	if(t) debTimer.reset();
}


void deb(int num){
	deb(num, false,'\n');
}


void deb(const char *format, ...){
	va_list args;
	va_start(args, format);
	_deb_(false, format, args);
	va_end(args);
}


void deb(int num, const char *format, ...){
	if(debugFg==0) return;
	deb(num,false,' ');
	va_list args;
	va_start(args, format);
	_deb_(false, format, args);
	va_end(args);
}
void debt(int num){
	deb(num,true,'\n');
}
void debt(const char *format, ...){
	va_list args;
	va_start(args, format);
	_deb_(true, format, args);
	va_end(args);
}


void debt(int num, const char *format, ...){
	if(debugFg==0) return;
	deb(num,false,' ');
	va_list args;
	va_start(args, format);
	_deb_(true, format, args);
	va_end(args);
}
void debt(){
	debTimer.reset();
}
//============================================================
//=======================      Timer   =======================
//============================================================
char *Timer::getTime(){
	long dt=getTimer();
	int ms=(int)(dt%1000);
	int s=(int)(dt/1000);
	snprintf(bb, sizeof(bb), "%i.%is",s,ms);
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


long mtime()
{
  struct timeval t;


  gettimeofday(&t, NULL);
  long mt = (long)t.tv_sec * 1000 + t.tv_usec / 1000;
  return mt;
}


char timerBufferQQ[256];
char *dateTime(){
	time_t lt=time(NULL);
	tm *t=localtime(&lt);
	snprintf(timerBufferQQ,sizeof(timerBufferQQ),"%02i.%02i.%02i %02i:%02i:%02i",t->tm_mday, t->tm_mon+1, t->tm_year%100,
			t->tm_hour, t->tm_min, t->tm_sec);
	return timerBufferQQ;
}


//============================================================
//====================   Files and Paths    ==================
//============================================================
int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}
//============================================================
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


//================ open file with control
FILE *xopen(const char* fname, const char *t){
	if(fname==0) errorExit("can\'t open file <null>");
	FILE* f=fopen(fname,t);
	if(f==0) errorExit("can\'t open file <%s>)",fname);
	return f;
}

// remove fucked backslash
char *correctFname(char* s){
	char *ss;
	for(ss=s;*ss;ss++) if(*ss=='\\') *ss='/';
	return s;
}
////================= create filename using path and name
//char* makeFileName(char *b, int siz, char *path, char*fname){
//	if(*fname=='~' || *fname=='/') return strcpy(b,fname);
//	if(*(lastChar(path))=='/')     snprintf(b, siz-1, "%s%s" ,path,fname);
//	else			               snprintf(b, siz-1, "%s/%s",path,fname);
//	return b;
//}
//================= create filename using path and name
char* makeFileName(char *b, char *path, char*fname){
	if(*fname=='~' || *fname=='/') path=0;
	if(path==0) return strcpy(b,fname);
	strcpy(b,path);
	if(*(lastChar(b)) != '/') strcat(b,"/");
	return strcat(b,fname);
//	if(*(lastChar(path))=='/')     sprintf(b, "%s%s" ,path,fname);
//	else			               sprintf(b, "%s/%s",path,fname);
	return b;
}
//================= create filename using path and name
char *makeFileName(char *b, char *path, char*fname, const char*ext){
	char bb[TBS];
	makeFileName(bb,path,fname);
	char *s=strrchr(bb,'/'); if(s==0) s=bb;
	char *sp=strrchr(s,'.'); if(sp  ) *sp=0;
	return strcat(strcat(b,"."),ext);
//	sprintf(b, "%s.%s",bb,ext);
//	snprintf(b, TBS, "%s.%s",bb,ext);
}
////================= create filename using path and name
//char *makeFileName(char *b, int siz, char *path, char*fname, const char*ext){
//	char bb[TBS];
//	makeFileName(bb,sizeof(bb), path,fname);
//	char *s=strrchr(bb,'/'); if(s==0) s=bb;
//	char *sp=strrchr(s,'.'); if(sp  ) *sp=0;
//	snprintf(b, siz-1, "%s.%s",bb,ext);
//	return b;
//}
//=================== extract fname wothout path
char *getFnameWithoutPath(char *buf, const char *fname){
	const char *s;
	s=strrchr(fname,'/'); if(s==0) s=fname; else s++;
	strcpy(buf,s);
	return buf;
}
//=================== extract fname wothout path and extension
char *getFnameWithoutExt(char *buf, const char *fname){
	getFnameWithoutPath(buf,fname);
	char *pp=strrchr(buf,'.');
	if(pp) *pp=0;
	return buf;
}
////================ Make Fname without path
//char *MakeFname(char *b, const char*fname, const char*ext){
//	char *s=getFnameWithoutExt(b,fname);
//	return strcat(strcat(s,"."), ext);
//}


//===================== platform independent Make Directory
int _makeDir(const char * path){
    struct stat sb;
    if (stat(path, &sb) == 0 && S_ISDIR(sb.st_mode)) return 0;
#if defined(_WIN32)
//#if __GNUC__ > 5
//	return _mkdir(path);
//#else
	return mkdir(path);
//#endif
#else
	mode_t mode=S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
	return mkdir(path, mode); // notice that 777 is different than 0777
#endif
}
//===================== platform independent Make Directory
void makeDir(const char *path){
	if(path==0 || *path==0) return;
	char b[TBS];
	strcpy(b,path);
	char *last=lastChar(b);
	if(*last=='/') *last=0;
	for(char *s=b; (s=strchr(s+1,'/'))!=0;){
		*s=0;
		if(_makeDir(b))
			errorExit("Can not create directory %s\n",b);
		*s='/';
	}
	_makeDir(b);
}
//====================== Lock file in platform - independent manner
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
	bool fg=false;				// The file do not exist. The header should be writen.
	FILE *f=fopen(fname,"r");	// check if file exists
	if(f!=0) {fg=true; fclose(f);}
	return fg;
}


////=================== Check if given file exists
//bool fileExists( char* path,  char *fname){
//	char b[TBS];
//	makeFileName(b, sizeof(b), path,fname);
//	return fileExists(b);
//}
////==================== check if the file exists
//bool fileExists( char* path,  char *fname, const char *ext){
//	char b[TBS];
//	makeFileName(b, sizeof(b), path,fname,ext);
//	return fileExists(b);
//}
//====================
unsigned long getFileTime(const char *fname){
	struct stat mystat;
	stat(fname,&mystat);
	return mystat.st_mtime;
}
//=================== extract file name
const char *getExt(const char *fname){
	const char *s=strrchr(fname,'/');
	if(s==0) s=fname;
	s=strrchr(s,'.');
	if(s==0) return 0;
	return s+1;
}



//============================================================
//========================    Memory    ======================
//============================================================
void zfree(void *a, const char* b){
	if(a) free(a); else if(b) writeLogErr("double free %s\n",b);
}
void *xmalloc(size_t n, const char *err){
	void *a=malloc(n);
	if(a==0){
		if(err==0)
			errorExit("can't allocate memory: %li",n);
		else
			errorExit("can't allocate memory. %s: %li",err,n);
	}
	return a;
}
void *xrealloc(void *a, size_t n, const char * err){
	a=realloc(a,n);
	if(a==0){
		if(err==0)
			errorExit("can't allocate memory: %li",n);
		else
			errorExit("can't allocate memory. %s: %li",err,n);
	}
	return a;
}


//============================================================
//========================    Random    ======================
//============================================================
//=====================================================
inline int longRand(){
	int x=rand();
	if(RAND_MAX > 0xfffff) {return x;}
	return (rand()<<15)|x;
}


double drand(){   /* uniform distribution, (0..1] */
	double x=(longRand()+1.0)/(LRAND_MAX+1.0);
  return x;
}
double drand(double xx){   /* uniform distribution, (0..1] */
	double x=(longRand())/(LRAND_MAX+1.0)*xx;
  return x;
}


int irand(int xx){   /* uniform distribution, (0..1] */
	long y=drand()*xx;
  return y;
}


double rGauss(){
	double phi=drand() * 2 * PI, r=0;
	while(r==0) r=drand();
	double rr=sqrt(-2.*log(r))*sin(phi);
	return rr;
}


double rExp(){
	return -log(drand());
}
// random gaussian variable with given mean and std deviation
double rGauss(double e, double sigma){
	return rGauss()*sigma+e;
}


// random integer in given interval
unsigned long randInt(unsigned long n){
	int k=longRand();
	double x=(double)(k)/(double)(LRAND_MAX);
	unsigned long rn=(unsigned long)(x*n);
	return rn;
}


//============================================================
//===================    Dynamic  Histogram    ===============
//============================================================
DinHistogram::DinHistogram(int ll){
	l=ll;
	getMem(hist[0],l,"Dinamic histogram 0");
	getMem(hist[1],l,"Dinamic histogram 1");
	getMem(cnts[0],l,"Dinamic histogram 3");
	getMem(cnts[1],l,"Dinamic histogram 4");
	clear();
}
DinHistogram::~DinHistogram(){
	xfree(hist[0],"Dinamic histogram 0");
	xfree(hist[1],"Dinamic histogram 1");
	xfree(cnts[0],"Dinamic histogram 3");
	xfree(cnts[1],"Dinamic histogram 4");
}


void DinHistogram::clear(){
	n[0]=n[1]=0;					//number of observations
	min=1.e+200; max=-min;			//min max values
	bin=hMin=hMax=0;				//bin width, max-min value in the histogram
	e[0]=e[1]=sd[0]=sd[1]=0;		//mean and std deviation
	zeroMem(hist[0],l);
	zeroMem(hist[1],l);
	zeroMem(cnts[0],l);
	zeroMem(cnts[1],l);
}




int DinHistogram::getIdx(double value){ //get index by value
	if(bin==0) return 0;
	if(value==hMax) return l-1;
	return (int)((value-hMin)/bin);
}
double DinHistogram::getValue(int idx){	//get value by index
	return hMin+idx*bin+bin/2;
}


double DinHistogram::getNormValue(int idx){
	double v=getValue(idx);
	v=(v-min)/(max-min);
	return v;
}




int DinHistogram::compress2Left(double value){  //Compress the histogram to the Left
	for(int i=0, j=0; i<l; i+=2, j++){			// Compress the values
		cnts[0][j]=cnts[0][i]+cnts[0][i+1];
		cnts[1][j]=cnts[1][i]+cnts[1][i+1];
	}
	for(int i=l/2; i<l; i++) cnts[0][i]=cnts[1][i]=0;	// Clear new space
	bin*=2;												// redefine bin size
	hMax=hMin+bin*l;									// redefine boundaries
	return getIdx(value);
}


int DinHistogram::compress2Right(double value){	//Compress the histogram to the Right
	for(int i=l-1,j=l-1; i > 0; i-=2,j--){		// Compress the values
		cnts[0][j]=cnts[0][i]+cnts[0][i-1];
		cnts[1][j]=cnts[1][i]+cnts[1][i-1];
	}
	for(int i=0; i<l/2; i++) cnts[0][i]=cnts[1][i]=0;	// Clear new space
	bin*=2;												// redefine bin size
	hMin=hMax-bin*l;									// redefine boundaries
	return getIdx(value);
}


void DinHistogram::addStat(double value, int count, int type){		//add the value to the statistics
	n[type]+=count; e[type]+=value; sd[type]+=value*value;	//count, mean, std dev.
	if(value > max) max=value;
	if(value < min) min=value;
}


void DinHistogram::add(double value, int type){			// add the value to the histogram
	add(value,1,type);
}


void DinHistogram::add(double value, int count, int type){			// add the value to the histogram
	if(count==0) return;
	if(n[0]==0 && n[1]==0){								// the histogram is empty
		cnts[type][0]=count; addStat(value,count,type);			// set the value
		hMin=hMax=value;								// define boundaries
		return;
	}
	if(hMin==hMax){										// the histogram contains only single bin
		if(value==hMin){								// new value equal to the single bin
			cnts[type][0]+=count; addStat(value,count,type);		//
			return;
		}
		else if(value > hMin){							// the new value is right to single bin
			cnts[type][l-1]+=count; addStat(value,count,type);		// the new value defines the right bin
			hMax=value;								// redefine the boundaries
		}
		else{											// the new value is less than single bin
			cnts[0][l-1]=cnts[0][0];					// Put the old value to the right end of the histogram
			cnts[1][l-1]=cnts[1][0];					//
			cnts[type][0]+=count; addStat(value,count,type);		// Put the new value to the left bin
			hMin=value;								// redefine the boundaries
		}
		bin=(hMax-hMin)/l;								// Define the bin
		return;
	}
	int i=getIdx(value);								// The histogram contains some values
	while(i < 0){i=compress2Right(value);}				// The new value is less than the first bin
	while(i >=l){i=compress2Left(value);}				// The new value is grater than the last bin
	cnts[type][i]+=count; addStat(value,count,type);				// Put the new value
}


void DinHistogram::fin(){								// Normalize and calculate the statistics
	for(int t=0; t<2; t++){
		for(int i=0; i<l; i++) hist[t][i]=(double) cnts[t][i]/n[t]/bin;
		e[t]/=n[t];
		sd[t]=(sd[t]-e[t]*e[t]*n[t])/(n[t]-1);
		sd[t]=sqrt(sd[t]);
	}
}


void DinHistogram::print(FILE* f){						// print the histogram
	fprintf(f,"#  min=%.3f max=%.3f \n",min,max);
	fprintf(f,"#  hMin=%.3f  hMax=%.3f bin=%.3f\n",hMin,hMax,bin);
	fprintf(f,"#  e0=%.3f sd0=%.3f n0=%i\n",e[0],sd[0],n[0]);
	fprintf(f,"#  e1=%.3f sd1=%.3f n1=%i\n",e[1],sd[1],n[1]);
	for(int i=getIdx(min); i<getIdx(max); i++){
		double h0=hist[0][i];
		double h1=hist[1][i];
		fprintf(f,"%.3f\t\%6i\t%.6f\t\%6i\t%.6f\n",getValue(i),cnts[0][i],h0,cnts[1][i],h1);
	}
}




//===================== Normalize the function to mean and sigma
double norm(double *x, int l){
	double d=0,e=0,dd,ee;
	for(int i=0; i<l; i++){d+=x[i]*x[i]; e+=x[i];}
	ee=e/l; d=d*l-e*e; dd=d/((l-1)*l);
	if(dd<0) {dd=0;} dd=sqrt(dd);
	if(dd <= ee*ee*1.e-5) {return 0;}
	for(int i=0; i<l; i++) x[i]=(x[i]-ee)/dd;
	return dd;
}


//======================================================================





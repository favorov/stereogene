/*
 * util.h
 *
 *  Created on: 05 дек. 2017 г.
 *      Author: andrey
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>

#ifndef UTIL_H_
#define UTIL_H_

const  int DEBUG_PRINT=1;
const  int DEBUG_LOG=2;

const  int BED_TRACK=1;
const  int BED_GRAPH=2;
const  int WIG_TRACK=3;
const  int BROAD_PEAK=4;
const  int MODEL_TRACK=8;
const  int CAGE_MIN=16;

const int LRAND_MAX=(RAND_MAX > 0XFFFFF)?RAND_MAX : 0x3fffffff;

extern unsigned long id;
extern char *logFileName;	// output filename
extern char *outFile;		// output filename
extern bool silent;				// inhibit stdout
extern bool verbose;				// number of suffle



const  double PI=3.14159265358979323846;
//================================= Timer ===================================
char *dateTime();
long mtime();
struct Timer{
	long start;
	char bb[80];
	Timer();
	void reset();
	long getTimer();
	char *getTime();
};

//=========================================================================
struct DinHistogram{		// Dynamic histogram for two variables
	int l;					//number of bins
	int    n[2];			//number of observations
	double min, max;		//min, max observed value
	double hMin, hMax;		// max-min value in the histogram
	double bin;				//bin width
	int    *cnts[2];		//counts
	double *hist[2];		//counts
	double e[2],sd[2];		//mean and std deviation

	DinHistogram(int ll);				//constructor
	~DinHistogram();
	int getIdx(double value);			// Get index for the value
	double getValue(int idx);			// Get value by index
	double getNormValue(int idx);		// Get normalized value (0..1)
	int compress2Left(double value);	// Compress the histogram to left
	int compress2Right(double value);	// Compress the histogram to left
	void addStat(double value, int count, int type);
	void add(double value, int type);				// Add the value
	void add(double value, int count, int type);				// Add the value
	void fin();
	void print(FILE* f);
	void clear();
};
//============================================ Memory operations
void *xmalloc(size_t n, const char * err);
void *xrealloc(void *a, size_t n, const char * err);
void zfree(void *a, const char* b);
//============================================ MACROS
#define getMem0(a,n,err) {if(a==0) a=(typeof a)xmalloc((n+100)*sizeof(*a),err);}
#define getMem(a,n,err)  {a=(typeof a)xmalloc((n+100)*sizeof(*a),err);}
#define del(a)  {delete a; a=0;}
#define realocMem(a,n,err)  {a=(typeof a)xrealloc(a,(n+100)*sizeof(*a),err);}
#define xfree(a,b) 		 {zfree(a,b); a=0;}
#define zeroMem(a,n) 	 {memset(a,0,n*sizeof(*a));}
#define xmemcpy(dest,pds, src, psrc, n) {memcpy(dest+(pds)*sizeof(*dest), src+(psrc)*sizeof(*src),(n)*sizeof(*src));}
#define max(a,b) a<b ? b : a
#define min(a,b) a>b ? b : a
#define abs(a) ((a) < 0 ? -(a) : (a))
#define sign(a) (a==0 ? 0 : (a<0 ? -1:1))
//=================================================== debug and error
extern const char *errStatus;
extern const char *debS;
extern int debugFg;


//=============================== Parsing ===========================
char *getAttr(char *s0, const char *name, char *buf); // find attribute by name in the string and return attribute value
char *skipSpace(char *s);		// skip space characters
char *skipNoSpace(char *s);		// skip non-space characters
char *strtoupper(char*s);		// convert string to upper-case
int  isEmpty(const char*buff);
bool isInt(const char *s);
bool isDouble(const char *s);
bool isUInt(const char *s);
char *trim(char *s);
int  keyCmp(const char *str, const char *key);
//=============================== Files =================================
FILE *xopen(const char*, const char*);		// open file if exists, exit otherwise
FILE *gopen(const char*, const char*);		// open file with parsing ~
bool fileExists(const char *fname);				// check if the file exists
bool fileExists(const char* path, const char *fname);				// check if the file exists
bool fileExists(const char* path, const char *fname, const char *ext); // check if the file exists
const char *getExt(const char *fname);					// extract file extension
char *getFnameWithoutExt(char *buf, char *fname);
void makeDir(const char *path);
unsigned long getFileTime(const char *fname);
char *correctFname(char* s);			// remove fucking MS Widows backslash
void flockFile(FILE *f);
void funlockFile(FILE *f);

//============================================ Arrays
double norm(double *x, int l);			// normalize to z-score

//============================================ Debugging and logging
void errorExit(const char *format, ...);
void clearLog();
void writeLog(const char *format, ...);
void writeLogErr(const char *format, ...);
void verb(const char *format, ...);
void xverb(const char *format, ...);
void deb(int num);
void deb(const char *format, ...);
void deb(int num, const char *format, ...);
void debt(int num);
void debt(const char *format, ...);
void debt(int num, const char *format, ...);
void debt(); 										//reset timer
void clearDeb();




#endif /* UTIL_H_ */

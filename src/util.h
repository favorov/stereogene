/*
 * util.h
 *
 *  Created on: 05 Dec. 2017
 *      Author: Mironov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>



#ifndef UTIL_H_
#define UTIL_H_

const int TBS=4096;			//(Temporary Buffer Size) size of buffers for filenames etc

const  int DEBUG_PRINT=1;
const  int DEBUG_LOG=2;


const  int BED_TRACK=1;
const  int BED_GRAPH=2;
const  int WIG_TRACK=3;
const  int BROAD_PEAK=4;
const  int NARROW_PEAK=5;
const  int MODEL_TRACK=8;
const  int CAGE_MIN=16;


const int LRAND_MAX=(RAND_MAX > 0XFFFFF)?RAND_MAX : 0x3fffffff;


extern unsigned long id;
extern char *logFileName;
extern char *repFile;		// reports filename
extern bool  silent;				// inhibit stdout
extern bool  verbose;				// number of suffle






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
struct Alias{
	char*	old;		// Pattern
	char*	rpl;		// Replace
	int 	l_rpl;		// rpl length
	int		l_old;		// old length

	Alias(const char *o, const char *r);
	int	replace(char *txt);
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
#define getMem0(a,n,err)  {if(a==0) a=(__typeof__(a))xmalloc((n+100)*sizeof(*a),err);}
#define getMem(a,n,err)   {a=(__typeof__(a))xmalloc((n+100)*sizeof(*a),err);}
#define getMemZ(a,n,err)  {a=(__typeof__(a))xmalloc((n+100)*sizeof(*a),err); memset(a,0,n*sizeof(*a));}
#define del(a)  {delete a; a=0;}
#define realocMem(a,n,err)  {a=(__typeof__(a))xrealloc(a,(n+100)*sizeof(*a),err);}
#define xfree(a,b) 		 {zfree(a,b); a=0;}
#define zeroMem(a,n) 	 {memset(a,0,n*sizeof(*a));}
//#define xmemcpy(dest,pds, src, psrc, n) {memcpy(dest+(pds)*sizeof(*dest), src+(psrc)*sizeof(*src),(n)*sizeof(*src));}
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
char * skipInt(char *s);
char *lastChar(char *s);
int  readInt(const char* val);
double readDouble(const char *s);
char *strtoupper(char*s);		// convert string to upper-case
int  isEmpty(const char*buff);
bool isInt(const char *s);
bool isDouble(const char *s);
bool isUInt(const char *s);
char *trim(char *s);
int  keyCmp(const char *str, const char *key);
bool isfloat(char *s);
//=============================== Files =================================
int isDirectory(const char *path) ;
FILE *xopen(const char*, const char*);		// open file if exists, exit otherwise
FILE *gopen(const char*, const char*);		// open file with parsing ~
bool fileExists(const char *fname);				// check if the file exists
//bool fileExists( char* path,  char *fname);				// check if the file exists
//bool fileExists( char* path,  char *fname, const char *ext); // check if the file exists
const char *getExt(const char *fname);					// extract file extension
void makeDir(const char *path);
unsigned long getFileTime(const char *fname);
char *correctFname(char* s);			// remove fucking MS Widows backslash
void flockFile(FILE *f);
void funlockFile(FILE *f);
//char *makeFileName(char *b, int siz, char *path, char*fname);	// make filename using path and name
//char *makeFileName(char *b, int siz, char *path, char*fname, const char*ext);	// make filename using path, name and extension
char *getFnameWithoutExt (char *buf, const char *fname);	// get fname without path & ext
char *getFnameWithoutPath(char *buf, const char *fname);    // get fname without path
char *makeFileName(char *b, char *path, char*fname);	// make filename using path and name
char *makeFileName(char *b, char *path, char*fname, const char*ext);	// make filename using path, name and extension
//char *MakeFname(char *b, const char*fname, const char*ext); // Make Fname without path and new ext


//============================================ Arrays
double norm(double *x, int l);			// normalize to z-score


//============================================ Debugging and logging
void errorExit(const char *format, ...);
void clearLog();
void writeLog(const char *format, ...);
void writeLogErr(const char *format, ...);
void verb(const char *format, ...);
void xverb(const char *format, ...);
//==================== just write debug information
void deb(int num);
void deb(const char *format, ...);
void deb(int num, const char *format, ...);
//===================== write debug info with time
void debt(int num);
void debt(const char *format, ...);
void debt(int num, const char *format, ...);
void debt(); 										//reset timer
void clearDeb();








#endif /* UTIL_H_ */

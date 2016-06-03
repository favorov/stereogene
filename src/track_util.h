/*
 * track_util.h
 *
 *  Created on: 17.02.2013
 *      Author: 1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>

#ifndef TRACK_UTIL_H_
#define TRACK_UTIL_H_

const  int DEBUG_PRINT=1;
const  int DEBUG_LOG=2;

#define PRM_EXT "prm"
#define BPROF_EXT "bprof"
#define AC_EXT ".acf"
#define DERIV '_'

#define  maxFactorCount        20

const  float NA=(1./1024./1024./1024./1024./1024.);
const  double PI=3.14159265358979323846;
extern const char* version;
const  int BED_TRACK=1;
const  int BED_GRAPH=2;
const  int WIG_TRACK=3;
const  int BROAD_PEAK=4;
const  int MODEL_TRACK=8;
const  int CAGE_MIN=16;

const  int BP_SCORE=0;
const  int BP_SIGNAL=1;
const  int BP_LOGPVAL=2;
const  int NONE=0;

const  int XML=1;
const  int TAB=2;

const  int GENE=4;
const  int IVS=0x10;
const  int EXON=0x20;
const  int GENE_BEG=0x41;
const  int GENE_END=0x42;
const  int IVS_BEG=0x11;
const  int IVS_END=0x12;
const  int EXON_BEG=0x21;
const  int EXON_END=0x22;

const  int KERN_NORM=1;
const  int KERN_LEFT_EXP =2;
const  int KERN_RIGHT_EXP=3;

const  int COLLINEAR=1;
const  int COMPLEMENT=2;
const  int IGNORE_STRAND=COLLINEAR|COMPLEMENT;
const  int MAX_SCORE=1;
const  int AV_SCORE=2;

const  int LOG_SCALE =1;
const  int LIN_SCALE =0;
const  int AUTO_SCALE=2;

const  int TOTAL=1;
const  int CHR_DETAIL=3;

const  int WIG_BASE=1;
const  int WIG_CENTER=2;
const  int WIG_SUM=0x10;
const  int WIG_MULT=0x20;


const int MAX_GENES=100000;

#define memSiz(a,n) n*sizeof(*a)
#define getMem0(a,n,err) if(a==0) a=(typeof a)xmalloc((n+100)*sizeof(*a),err)
#define getMem(a,n,err) a=(typeof a)xmalloc((n+100)*sizeof(*a),err)
#define zeroMem(a,n) memset(a,0,n*sizeof(*a))
#define xfree(a) {if(a) free(a); a=0;}
#define max(a,b) a<b ? b : a
#define min(a,b) a>b ? b : a
#define abs(a) (a<0 ? -a : a)

struct Model;
//========================= Common parameters ================================
extern int  stepSize;          // frame size of profile
extern long long GenomeLength;      // TOTAL LENGTH OF THE GENOME
extern char *chromFile;		   // Chromosomes file name
extern char *profile1;		   // first profile file file name
extern char *profile2;		   // second profile file file name
extern char *trackFil;		   // input track
extern char *aliaseFil;		   // aliases
extern char *mapFil;		   // map file
extern char *inputProfiles;    // input formula for track combination
extern int 	bpType;

extern bool  writeDistr;		// flag: write distributions

extern char *trackName;   	// current track name
extern char *profPath;		// path to binary profiles
extern char *trackPath;		// path to GBrowse track files
extern char *resPath;		// path to results files (tracks and distributions)
extern char *outFile;		// output filename
extern char *logFileName;	// output filename
extern char *statFileName;	// File name for cummulative statistics
extern char *paramsFileName; // Filename for save parameters of runs
extern unsigned long id;

//============================ Track parameters ===============================
extern int   complFg;		// Flag: IGNORE_STRAND - ignore strand; COLLINEAR - compare collinear strands; COMPLEMENT - compare complement strands
extern int   intervFlag0;	// Flag: for bed tracks take intervals with score 1.
extern bool  NAFlag;		// Flag: 0 -> uncovered_values=0; otherwise uncovered values=NA
extern bool  strandFg0;	 	// Flag: 1 - strand-dependent
extern int   profileLength;	// size of the profile array
extern float *profile;		// uncompressed profile array
extern float *profilec;		// uncompressed profile array

extern unsigned char *byteprofile;	 // compressed profile array
extern unsigned char *byteprofilec;	 // compressed profile array for complement strand
extern int   logScale;		// use logScale
extern double scaleFactor0;

extern char *pcorProfile;    	// partial correlation profile file name
extern float *outTrackProfile;  // correlation track

//=========================== Correlation parameters ===========================

extern int flankSize;    // size of flank (nucleotides)
extern int LFlankProfSize;        // size of the left flank (profile scale)
extern int RFlankProfSize;        // size of the right flank (profile scale)
extern double noiseLevel; // level of noise
extern int wSize;        // size of widow (nucleotides)
extern int wStep;        // window step   (nucleotides)
extern int wProfSize;    // size of widow (profile scale)
extern int wProfStep;    // window step (profile scale)
extern double kernelSigma;    	// kernel width (nucleotides)
extern double kernelShift;    	// Kernel mean (for Gauss) or Kernel start for exponent
extern double kernelNS;			// Correction for non-specifisity
extern double kernelProfSigma;  // kernel width (profile scale)
extern double kernelProfShift;  // Kernel shift (profile scale)


extern int kernelType;   		// kernel type
extern double maxNA0;				// max allowed NA values - if at lest one of profiles contains more than maxNA NA values the window will be ignored
extern double maxZero0;				// max allowed zero values - if both profiles contains more than maxZero zeros the window will be ignored
extern double maxNA;				// max allowed NA values - if at lest one of profiles contains more than maxNA NA values the window will be ignored
extern double maxZero;				// max allowed zero values - if both profiles contains more than maxZero zeros the window will be ignored
extern int nShuffle;			// number of shuffle in percents
extern int maxShuffle;		    // max number of shuffle
extern int minShuffle;		    // min number of shuffle
extern int nCompare;			// number of observations

extern double pVal;				// min (-log10(p-value)) for output
extern double qVal;				// min (-log10(q-value)) for output
extern bool verbose;				// number of suffle
extern int threshold;
extern int lAuto;
extern int lProfAuto;
extern int corrScale;			// scale for correlations
extern bool corrOnly;			// only calculate corr functions
extern bool writeBPeak;			// write BroadPeak
extern bool silent;				// inhibit stdout
extern bool syntax;				// Strong syntax control

extern int writeDistCorr;		// write BroadPeak

extern int outWIG;

extern double prod11,prod12,prod22,sprod11, sprod12,sprod22;	//cummulative scalar products
extern int profWithFlanksLength; // size of profWindow array (including random flanks)
extern int nProd;
extern int pcaFg;
extern int nPca;
extern int pcaSegment;	//segment length in profile scale
extern double totCorr;
extern bool RScriptFg;
extern bool outSpectr;
extern int outRes;
extern bool outChrom;
extern int outThreshold;
extern const char *errStatus;
extern const char *debS;
extern int debugFg;
extern double *xDat,*yDat,*xyCorr;  	//Woriking arrays

extern int cage;
extern bool clearProfile; //Force profile recalculation

extern int scoreType;

extern int nBkg, nFg;					// size of background and foreground sets of data
extern double *BkgSet, *FgSet;			// background and foreground sets of the correlations


extern FILE *logFile;


void deb(int num);
void deb(const char *format, ...);
void deb(int num, const char *format, ...);
void clearDeb();


struct Chromosome{
    char *chrom;		//Chromosome name
    long length;		//Chromosome length
    int  base;			//Start position in binary profile

    double av1, av2, corr, lCorr, count;
    int densCount;
    float *distDens;
    Chromosome(char *chr, long l, int bb);
    void clear();
};

struct GenomePos{
    char* chrom;	// Chromosome name
    long pos;		// Chromosom position
};

struct ScoredRange{
	Chromosome *chr;
    char* chrom;	// Chromosome name
    long beg;		// start Chromosome position
    long end;		// end chromosome position
    float score;	// score
    ScoredRange();
};
//              mapSgm=b:-500..e:+200
struct FileName{
	const char *path;	// Path
	const char *name;	// File name
	const char *ext;	// extension
	FileName(char *fn);	// constructor parses name
	char *fname();		// generate full filename
	char *coreFname();  // fileneme without extension
};
//===============================================================
struct MapRange{
	int f,t;
	int cumLength;
	MapRange();
	MapRange(int f, int t);
};

struct IVSet{
	MapRange **ivs;
	int capacity;
	int nIv;
	int totLength;
	int ivNo;
	IVSet();
	void addIv(int f, int t);
	int randPos();
	void fin();
	void clear();
	void write(FILE*f);
	void print(int f, int t);
};

//===================================================================
class Exon{
public:
	long from;
	long to;
};

class Exons{
public:
	Exon *list;
	int n;
};

//===============================================================
struct bTrack{		        // Binary track
	unsigned char *bytes;	// byte_track
	unsigned char *cbytes;	// complement byte_track
	double *profWindow;		// decoded profile with flanks
	char *name;				// track name
	int lProf,		// profile length (size of bytes array), step (in profile scale)
		lScale;		// logScale (0 = linear scale, 1 = log scale)
	double  av,				// average score
			sd, 			// score standard deviation
			minP, 			// min score value
			maxP, 			// max score value
			bScale, 		// scale factor
			delta,			// =min-av
			av0,
			sd0,
			nn;

	bool strandFg;
	float scaleFactor;

	int intervFlag;
	bool hasCompl;
	IVSet ivs, ivsC;
	int trackType;
	int deriv;						// order of derivative
	double projCoeff;

	bTrack();						// empty constructor
	bTrack(const char* fname);		// constructor that reads file
	void read(const char *fname);	// read file
	double *getProfile(int pos, bool cmpl);	// decode profile starting with given position
	double getVal(unsigned char b);
	double getValue(int pos, int cmpl);
	double getProjValue(int pos, bool cmpl);
	double addStatistics();
	void finStatistics();
	int readProfileToArray(double *x, int scale, int from, int to, bool cmpl);
	int  countNA  (int pos, bool cmpl);		// count NA elements in the window
	int  countZero(int pos, bool cmpl);		// count zero elements in the window

	void printBytes(int from, int to); // Print bytes (for testing)
	void printWindow(int n);
	void printBytes(FILE*f,int from, int to); // Print bytes (for testing)
	void writeBytes(FILE*f); // Print bytes (for testing)
	void printWindow(FILE*f,int n);
	void correct(bTrack bPC);
	void clear();
	void ortProject();				// calculate projection coeff
	void makeIntervals();
	void makeMapIntervals();
	void makeIntervals(unsigned char *bytes, IVSet *iv);
	void makeMapIntervals(unsigned char *bytes, IVSet *iv);
	int getRnd(bool cmpl1);
	bool check(const char *fname);
	void makeBinTrack();
	void initProfile();
	void readTrack(const char *fname, int cage=0); //cage: if cage>0 : end=beg+cage; cage<0: beg=end+cage.
	void trackDef(char *s);
	void trackAutoCorrelation();
	void finProfile();
	float normProf(float x);
	void writeByteProfile();
	void writeProfilePrm();
	void addSgm(ScoredRange *bed, float *prof);
	void addSgm(char strnd, ScoredRange *bed, float *prof, float *profc, int strndFg);
};
//==============================
struct Term{
	char * fname;
	float mult;
	int read(char *b);
	void add();
	void make();
};

struct Model:bTrack{
	char * definition;
	Term trm[120];
	int nTerm;

	Model();
	void readMap(char *fnam);
	void create();
	void write();
};
struct Timer{
	long start;
	char bb[80];
	long mtime();
	Timer();
	void reset();
	long getTimer();
	char *getTime();
};

#define  maxPrimeFactor        11
#define  maxPrimeFactor2       maxPrimeFactor*2
#define  maxPrimeFactorDiv2    (maxPrimeFactor2+1)/2
#define  maxFactorCount        20

//===============================================================
class Fourier{		// Fourier transformation
public:
	int err;
	int length;				// array length
	double *datRe,*datIm, 			// im part of input data (always 0)
			*re, *im; 		// real and im parts of transformation
	double re0,im0;

	Fourier();				// empty constructor
	Fourier(int n);			// constructor
	~Fourier();
	void freeMem();
	void init(int len);
	void setDat(double *reD);
	void setDat(double *reD, double *imD);
	void calc(double *reD, double *imD, int deriv);	//do transform with real and Im data
	void calc(double *dat, int deriv);     // Do transformation with given real data
	void calc0(double *dRe, int deriv);	   // Do transform and set re[0]=im[0]=0
	void calc(int deriv); 			       // Do transformation
	void derivat();						   // get derivative from fft
	void restore(){re[0]=re0; im[0]=im0;}
};

struct Complex{
	double re,im;
	Complex(){re=0; im=0;}
	double Mod();
	Complex scalar(Complex otherC);
};

struct statTest{		    // result of a statistical test
	int n1,n2;				// sample sizes
	double u,e,sigma,z,pVal;	// statistics, mean, deviation, z-score and p-value for statistics
};

struct Kernel{				// Generic kernel
public:
	char *name;				// kernel name
	int length;				// length (equal to the profile length
	double * kern;			// the kernel values
	double * ckern;			// the complement kernel values
	Fourier ft;				// Fourier transformation for the kernel
	Fourier cft;			// Fourier transformation for the complement kernel
	Fourier fx,fy;			// Fourier transformation for the input data
//	Fourier fpc;			// Fourier transformation for partial correlation track
	bool hasCompl;			// for symmetrical Kernel flag=false; otherwise flag=1;
	virtual ~Kernel(){;}
	void init(int l);
	void fft();
	void fftx(double* , int deriv);		// do transform for given data
	void ffty(double* , int deriv);		// do transform for given data
	void fftpc(double*);		// do transform for given data
	double scalar(Fourier *f1, Fourier *f2, Complex *c, bool complem);
	double dist  (Fourier *f1, Fourier *f2, bool complem);
	double dist  (Fourier *f1, Fourier *f2, Fourier *fpc, bool complem);
	double dist(bool complem);
	void makeKernel(int l);
	double NSCorrection(double x, double val);
	virtual double kernVal(double x)=0;
	void restore(){fx.restore(); fy.restore();}
};

struct  LeftExpKernel:Kernel{
public:
	double sigma, e;
	LeftExpKernel();
	LeftExpKernel(double e,double sgm, int l); //Nucleotides
	double kernVal(double x);
};

struct  RightExpKernel:Kernel{
public:
	double sigma, e;
	RightExpKernel();
	RightExpKernel(double e,double sgm, int l); //Nucleotides
	double kernVal(double x);
};


struct  NormKernel:Kernel{
public:
	double sigma, e;
	NormKernel();
	NormKernel(double e,double sgm, int l);	//Nucleotides
	double kernVal(double x);
};

struct Histogram{
	double  minVal, maxVal,  // Min & Max values. Min=-1; Max=1;
			bin,			 // bin size
			e,  			 // Mean
			sigma,			 // standard deviation
			beta;			 // parameter fof Beta-distribution
	int 	nBin,			 // Number of bins
			count;			 // Number of observations
	double  *dd,			 // Distribution density
			*db,			 // Appropriate Beta density
			*Fp,			 // Cummulative  left distribution  ( F(x) )
			*Fm;			 // Cummulative  right distribution ( 1-F(x) )
	int iq,im,iqq;

	bool ready;				 // block add() and norm() after norm()
	Histogram(int nBin);
	void add(double x);
	void norm();
	void normBeta();		 // Calculate cummulative Beta distribution
	void normF();			 // Calculate cummulative real distribution
	double pValp(double x);
	double pValm(double x);
	double interpol(double x, double* fun);
	void print(FILE *f);
	void fitBeta();
	double error(double b);
};

struct MapPos{
	bool fg;				// flag: true-> position relative to the begin (upstream)
	int pos;				// position
	bool scaled;
	MapPos();
	int read(char *s);
	void scale(int pstep);
	void print();
	void print(char *b);
};

struct Interval{
	int f,t;
};

struct MapIv{				// description of the map interval
	MapPos beg,end;

	MapIv();
	void scale(int pstep);
	int read(char* s);		// read interval. Syntax: b:-500..b:200 = upstream region
	void print();
	char* print(char *b);
};

struct Aliase{
	char *oldName;
	char *newName;
	int lnew, lold;
};
class AliaseTable{
public:
	Aliase *als;
	int nAls;
	AliaseTable(){nAls=0; als=0;}
	void readTable(const char* fname);
	char *convert(char*oldName);
};

struct FileListEntry{
	int id;
	char *fname;
};

struct Correlation{
	float *correlation, *corrPlus, *corrMinus, *spectrumX, *spectrumY;  	//correlation function
	int nCorr, nPlus, nMinus;				//counter
	double min,max,av,sd;
	Correlation();
	void init();
	void getLimits(int &left, int &right, double &bottom, double &top);
	void norm();
	void calcWindowCorrelation(int pos, bool cmpl1, bool cmpl2,double corr);
//	void print(char *fname);
	void printSpect(char *fname);
};

struct PairEntry{			//==== correlation for pair of windows
	int profPos;			//==== Profile position
	float d;				//==== Correlation
};


//=====================================================================
extern AliaseTable alTable;
extern Chromosome *chrom_list;       // list of chromosomes
extern int n_chrom;
extern Chromosome *curChrom;
extern bTrack bTrack1, bTrack2, projTrack, mapTrack; // Binary profiles
extern Kernel *kern;			// current kernel
extern MapIv miv;				// interval for mapping
extern Model model;
extern FileListEntry files[256];
extern int   nfiles;

extern Correlation correlation;		    // array for correlation picture
extern Correlation bgcorrelation;		// array for background correlation picture

extern Fourier wCorrelation;			// distance correlation calculation

extern PairEntry *pairs;				// array for pair's correlation (foreground)
extern int nPairs;						// number of foreground observations
extern Histogram bgHist;				// Background Histogram initiation
extern Histogram fgHist;				// Foreground Histogram initiation
extern double BgAvCorr;					// average Background correlation
extern double FgAvCorr;					// average Foreground correlation

//=============================== Chromosomes ===========================
int  readChromSizes(char *fname);			// read chromosomes
Chromosome *getChromByPos(int pos);			// get Chromosome by file position
long pos2filePos(char*chrom,long pos);		// transform genome position to byte-profile position
void filePos2Pos(int pos, ScoredRange *gr, int length); // transform byte-profile position to genome position
Chromosome* findChrom(char *ch);			// Find chromosome record by name
Chromosome *checkRange(ScoredRange *gr);    // check if genome range is valid
void clearChromosomes();
void addChromStat(int pos, bool cmpl1,bool cmpl2,double corr);
//=============================== Parsing ===========================
char *getAttr(char *s0, const char *name, char *buf); // find attribute by name in the string and return attribute value
char *skipSpace(char *s);		// skip space characters
char *skipNoSpace(char *s);		// skip non-space characters
char *strtoupper(char*s);		// convert string to upper-case
int EmptyString(const char*buff);
bool isInt(const char *s);
bool isDouble(const char *s);
bool isUInt(const char *s);
//=============================== File names ===========================
void makeDir(const char *path);
char *makeFileName(char *b, const char *path, const char*fname);	// make filename using path and name
char *makeFileName(char *b, const char *path, const char*fname, const char*ext);	// make filename using path, name and extension
char *correctFname(const char* s);			// remove fucking MS Widows backslash
char *cfgName(char* p, char* ext);			// Make config file name
char *makePath(char* pt);					// Make path - add '/' to the end of pathname
FILE *xopen(const char*, const char*);		// open file if exists, exit otherwise
FILE *gopen(const char*, const char*);		// open file with parsing ~
char *getFname(char *s);					// get filename without path
bool fileExists(const char *fname);				// check if the file exists
bool fileExists(const char* path, const char *fname);				// check if the file exists
bool fileExists(const char* path, const char *fname, const char *ext); // check if the file exists
void makeDirs();
const char *getExt(const char *fname);					// extract file extension
char *getFnameWithoutExt(char *buf, char *fname);
int   getTrackType(const char *fname);
unsigned long getFileTime(const char *fname);
void flockFile(FILE *f);
void funlockFile(FILE *f);
void addFile(const char* fname);



//============================================== read config file
char *trim(char *s);
//void readArg(char *b);
//int readArg(char *s1, char *s2, int check);
//void readArgs(int argc, const char *argv[]);
//void readCfg(int argc, const char *argv[]);
//void readCfg(char *cfg);
int  keyCmp(const char *str, const char *key);
int  getFlag(char*s);
const char*getKernelType();
const char*getPC();//partial correlation
void  makeId();
const char *getIvFlag();

//============================================== random & statistics
double rGauss();							// standard normal random
double rGauss(double e, double sigma);		// normal random with given mean and sigma
unsigned long randInt(unsigned long n);							// uniform random int
double xBetaD(double beta1, double beta2, double x);	// beta density
statTest *MannWhitney( double *set1, int nSet1, double *set2,int nSet2);	//Mann-Whitney test

//============================================ Fourier transformation
int nearPow2(int n, int &i);
int nearPow2(int n);
int nearFactor(int n);
void printFactors(int l);
extern "C" void fftl(int n,double are[],double aim[],double bre[],double bim[]);
extern "C" void fft(double xRe[], double xIm[],double yRe[], double yIm[]);
extern "C" void initFFT(int n);

//============================================ Correlations
int corrFunc(double *x, double *y, double *rc, int l); //calculate autocorrelation and correlation
				//====== x - first variable (x autocorr. will be stored here),
				//====== y - second variable (y autocorr. will be stored here),
				//====== rc - array for correlation,
double storeCorrTrack(int pos, bool cmpl1, bool cmpl2);
void initOutWig();
void finOutWig();
void printStat();
void printFgDistr();
void printBroadPeak();
void printR();
void printCorrelations();
void printChomosomes(char *fname);
void printChrDistances(char *fname);
int  Correlator();
int  Preparator(const char *fname);

//============================================ Arrays
double arrayMax(double *d, int n);		// find max value
double norm(double *x, int l);			// normalize to z-score
void *xmalloc(size_t n, const char * err);
void errorExit(const char *format, ...);
void writeLog(const char *format, ...);
void verb(const char *format, ...);
void xverb(const char *format, ...);
void helpPage();
//======================== Mapping
void mapIntervals();
void fillMap();

//======================== Testing procedures
void printBytes(unsigned char* byteprofile,int from, int to);
void printProfile(float *profile, int from, int to);
void printProfile(FILE *f, float *profile, int from, int to);
void printChrom();
void testFourier();
void testFourierPC();
void GenerateData(FILE *bed, FILE *wig, Chromosome *ch);
void GenerateData();
int genmain();
double *readDistr(const char *fname, int &nn);
void test();
void pcaMain(const char *fname);
void CageMin(const char *fname1, const char *fname2);
void Active();
#endif /* TRACK_UTIL_H_ */

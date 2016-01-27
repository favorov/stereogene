#include "track_util.h"
#include <unistd.h>

#if defined(_WIN32)
#include <conio.h>
#else
#include <termios.h>
#endif


int xpause(){
	int c=0;
#if defined(_WIN32)
	c=_getch(); printf("\n%c\n",c);
#else
	return 0;
	struct termios oldt, newt;
	tcgetattr( STDIN_FILENO, &oldt);
	newt = oldt;
	newt.c_lflag &= ~(ICANON);
	tcsetattr( STDIN_FILENO, TCSANOW, &newt);
	c=getchar();
	tcsetattr( STDIN_FILENO, TCSANOW, &oldt);
#endif
	return c;
}

const int PRM_INT=1;
const int PRM_DOUBLE=2;
const int PRM_STRING=3;
const int PRM_ENUM=4;
const int PRM_FG=5;
const int PRM_MAPIV=6;

const int PRM_UNKNOWN=-0XFFFFFFF;

char * cfgFile=0;

struct Name_Value{			// symbolic name for a value
	const char* name;		// name for the value
	int value;				// value
	Name_Value(const char *nm, int val){name=nm; value=val;}
};

struct Param{
	const char* name;			// command line (cfg) argument name
	int type;					// parameter type
	Name_Value **enums;			// array of aviable values
	void *prm;					// pointer to the argument value
	int value;					//a value that should be set if -prm is used in a command line
	const char *description;
	Param(const char* descr);
	Param(const char* _name, int    *prm, const char* descr);
	Param(const char* _name, int    *prm, int val, const char* descr);
	Param(const char* _name, int    *prm, Name_Value **fg, const char* descr);
	Param(const char* _name, bool   *prm, const char* descr);
	Param(const char* _name, bool   *prm, bool val, const char* descr);
	Param(const char* _name, double *prm, const char* descr);
	Param(const char* _name, char * *prm, const char* descr);
	Param(const char* _name, MapIv  *prm, const char* descr);
	void setVal();
	void init(const char* _name, void* _prm, int type, Name_Value **fg, const char* descr);
	void printDescr();
	int readVal(char *s);
	int readEnum(char *s);
};

//===================================================================================================
Param *findParam(const char * name);
Param * readParam(char *name, char *value);
Param *findParam(const char * name);
void parseArgs(int argc, char **argv);
void readPrm(char *s);
void readPrm(char *key, char *val);
void printHelp();


//==============================================================================
Name_Value* intervalTypes[]={
		new Name_Value("NONE",NONE),
		new Name_Value("GENE",GENE),
		new Name_Value("EXON",EXON),
		new Name_Value("IVS",IVS),
		new Name_Value("GENE_BEG",GENE_BEG),
		new Name_Value("EXON_BEG",EXON_BEG),
		new Name_Value("IVS_BEG",IVS_BEG),
		new Name_Value("GENE_END",GENE_END),
		new Name_Value("EXON_END",EXON_END),
		new Name_Value("IVS_END",IVS_END),
		0
};
Name_Value* scaleTypes[]={
		new Name_Value("LOG",LOG_SCALE),
		new Name_Value("LIN" ,LIN_SCALE),
		new Name_Value("AUTO",AUTO_SCALE),
		0
};
Name_Value* bpTypes[]= {
		new Name_Value("SCORE",BP_SCORE),
		new Name_Value("SIGNAL" ,BP_SIGNAL),
		new Name_Value("LOGPVAL",BP_LOGPVAL),
		0
};
Name_Value* kernelTypes[]= {
		new Name_Value("NORMAL",KERN_NORM),
		new Name_Value("LEFT_EXP",KERN_LEFT_EXP),
		new Name_Value("RIGHT_EXP",KERN_RIGHT_EXP),
		0
};

Name_Value* complFlags[]={
		new Name_Value("IGNORE_STRAND",IGNORE_STRAND),
		new Name_Value("COLLINEAR",COLLINEAR),
		new Name_Value("COMPLEMENT",COMPLEMENT),
		0
};
Name_Value* outWigTypes[]={
		//	correlation is ( f*\int g\rho + g*\int f\rho )
		new Name_Value("NONE",NONE),
		new Name_Value("BASE",WIG_BASE|WIG_SUM), //correlation without substract average
		new Name_Value("CENTER",WIG_CENTER|WIG_SUM),//substract average
		//	correlation is ( \int g\rho * \int f\rho )
		new Name_Value("BASE_MULT",WIG_BASE|WIG_MULT),
		new Name_Value("CENTER_MULT",WIG_CENTER|WIG_MULT),
		0
};
Name_Value* distTypes[]={
		new Name_Value("TOTAL",TOTAL),
		new Name_Value("DETAIL",CHR_DETAIL),
		new Name_Value("NONE",NONE),
		0
};
Name_Value* outResTypes[]={
		new Name_Value("NONE",NONE),
		new Name_Value("XML",XML),
		new Name_Value("TAB",TAB),
		new Name_Value("BOTH",XML|TAB),
		0
};

//===================================================================================================
Param *pparams[]={
//================================================== Common parameters
		new Param("common parameters"),
		new Param("v"		    , &verbose		,1, "verbose"),
		new Param("verbose"	    , &verbose		,"verbose"),
		new Param("s"		    , &silent		,1, "no output to stdout"),
		new Param("silent"	    , &silent		,"no output to stdout"),


//=============================================================================================================
		new Param("preparation parameters"),
		new Param("step" 	 	, &stepSize  	,"step for input averaging"),
		new Param("scale" 	 	, &logScale  	,scaleTypes,"use Logarithmic or linear scale"),
		new Param("scaleFactor"	, &scaleFactor0 ,0),
		new Param("clear" 	  	, &clearProfile ,"force binary profile preparation"),
		new Param("c" 	  	    , &clearProfile , 1,"force  binary profile preparation"),

//=============================================================================================================
		new Param("paths and files"),
		new Param("cfg" 		, &cfgFile 		,"config file"),
		new Param("profPath" 	, &profPath 	,"path for binary profiles"),
		new Param("trackPath" 	, &trackPath 	,"path for tracks"),
		new Param("resPath" 	, &resPath 		,"path for results"),

		new Param("statistics" 	, &statFileName	 ,"cumulative file with statistics"),
		new Param("params" 		, &paramsFileName,"cumulative file with parameters"),
		new Param("log" 		, &logFileName	 ,"cumulative log-file"),
		new Param("aliases" 	, &aliaseFil	,0),
//=============================================================================================================
		new Param("input parameters"),
		new Param("chrom"		, &chromFile	,"chromosome file"),
		new Param("strand" 	 	, &strandFg0  	,"account for strand information"),
		new Param("intervals"	, &intervFlag0	,intervalTypes,"interval type in BED file"),
		new Param("gene"		, &intervFlag0	,GENE		,"consider entire gene"),
		new Param("exon"		, &intervFlag0	,EXON		,"consider exons"),
		new Param("ivs"			, &intervFlag0	,IVS 	 	,"consider introns"),
		new Param("gene_beg"	, &intervFlag0	,GENE_BEG	,"gene starts"),
		new Param("exon_beg"	, &intervFlag0	,EXON_BEG	,"exons starts"),
		new Param("ivs_beg"		, &intervFlag0	,IVS_BEG	,"introns starts"),
		new Param("gene_end"	, &intervFlag0	,GENE_END	,"gene ends"),
		new Param("exon_end"	, &intervFlag0	,EXON_END	,"exons ends"),
		new Param("ivs_end"		, &intervFlag0	,IVS_END	,"introns ends"),
		new Param("bpType" 	  	, &bpType  		,bpTypes	,"The value used as a score for BroadPeak input file"),
		new Param("pcorProfile" , &pcorProfile	,"Track for partial correlation"),
		new Param("NA"       	, &NAFlag     	,1 , "use NA values as unknown and fill them by noise"),
		new Param("threshold"	, &threshold	,"threshold for input data for removing too small values"),
		new Param("map" 		, &mapFil 		,0),
		new Param("mapIv" 		, &miv			,0),
//=============================================================================================================
		new Param("Analysis parameters"),
		new Param("kernelType"	, &kernelType	,kernelTypes,0),
		new Param("kernelSigma"	, &kernelSigma  ,"Kernel width"),
		new Param("kernelShift"	, &kernelShift  ,0),
		new Param("wSize" 	  	, &wSize  		,"Window size"),
		new Param("wStep" 	  	, &wStep  		,0),
		new Param("kernelNS"	, &kernelNS  	,0),
		new Param("flankSize"	, &flankSize  	,0),
		new Param("maxNA"		, &maxNA  		,"Max number of NA values in window (percent)"),
		new Param("maxZero"		, &maxZero  	,"Max number of zero values in window (percent)"),
		new Param("nShuffle"	, &nShuffle  	,"Number of shuffles for background calculation (percent of window pairs)"),
		new Param("MaxShuffle"	, &maxShuffle  	,"Max number of shuffles"),
		new Param("MinShuffle"	, &minShuffle  	,"Min number of shuffles"),
		new Param("noiseLevel"	, &noiseLevel  	,0),
		new Param("complFg"		, &complFg		,complFlags,0),
//=============================================================================================================
		new Param("Output parameters"),
		new Param("outBPeak" 	, &writeBPeak   ,"write BradPeak file"),
		new Param("pVal"		, &pVal  		,"threshold for BroadPeak output"),
		new Param("qVal"		, &qVal  		,"threshold for BroadPeak output"),
		new Param("outSpectr" 	, &outSpectr   	,"write fourier spectrums"),
		new Param("outChrom" 	, &outChrom   	,"write statistics by chromosomes"),
		new Param("writeDistr" 	, &writeDistr   ,"write foreground and background distributions"),
		new Param("Rscrpit" 	, &RScriptFg   	,0),
		new Param("r" 			, &RScriptFg   	,1,"write R script for the result presentation"),
		new Param("Distances" 	, &writeDistCorr,distTypes	,"Write distance correlations"),
		new Param("outWig"		, &outWIG		,outWigTypes,"parameters for local correlation file"),
		new Param("outThreshold", &outThreshold	,"threshold for output to correlation profile"),
		new Param("corrOnly" 	, &corrOnly   	,0),
		new Param("corr" 		, &corrOnly   	,1, 0),
		new Param("lAuto"	  	, &lAuto  		,0),
		new Param("outRes" 		, &outRes 		,outResTypes,"format for results in statistics file"),
		new Param("pcaSegment" 	, &pcaSegment	,0),
		new Param("nPca" 		, &nPca			,0),
		new Param("cage" 		, &cage			,0),

		new Param("Happy correlations!"),
		0,
};

//=================================================================================================
//===================================  End declaration ============================================
//=================================================================================================
void Param::init(const char* _name, void *_prm, int _type, Name_Value **fg, const char* descr){
	name=_name;
	enums=fg;
	description=descr;
	prm=_prm;
	type=_type;
	value=PRM_UNKNOWN;
}
Param::Param(const char* descr)												{init("",0,0,0 ,descr);}
Param::Param(const char* _name,  int    *_prm, const char* descr)			{init(_name,_prm,PRM_INT	,0 ,descr);}
Param::Param(const char* _name,  int    *_prm, Name_Value **fg, const char* descr){init(_name,_prm,PRM_ENUM	,fg,descr);}
Param::Param(const char* _name,  bool   *_prm,  const char* descr)			{init(_name,_prm,PRM_FG		,0 ,descr);}
Param::Param(const char* _name,  double *_prm,  const char* descr)			{init(_name,_prm,PRM_DOUBLE	,0 ,descr);}
Param::Param(const char* _name,  char*  *_prm, const char* descr)			{init(_name,_prm,PRM_STRING	,0 ,descr);}
Param::Param(const char* _name, MapIv   *_prm, const char* descr)			{init(_name,_prm,PRM_MAPIV	,0 ,descr);}
Param::Param(const char* _name,  int    *_prm, int val, const char* descr)	{init(_name,_prm,PRM_INT	,0 ,descr);
	value=val;}
Param::Param(const char* _name,  bool    *_prm, bool val, const char* descr)	{init(_name,_prm,PRM_FG	,0 ,descr);
	value=val;}
//====================================================================================================
int Param::readEnum(char *s){
	for(int i=0; enums[i]!=0; i++){
		if(keyCmp(enums[i]->name, s) ==0) return enums[i]->value;
	}
	return PRM_UNKNOWN;
}

//============================ return 0 -> OK ==========================================================
int Param::readVal(char *s){
	s=trim(s);
	int fg=0;
	switch(type){
	case PRM_INT:  		{
		s=strtok(skipSpace(s),".,; \t");
		if(isInt(s))
			*(int *)  (prm)=atoi(s);
		else
			fg=1;
		break;
	}
	case PRM_FG:		{
		int x=getFlag(s);
		if(x==-1) fg=1;
		else *(bool*)  (prm)=x;
		break;
	}
	case PRM_DOUBLE:   	{
		s=strtok(skipSpace(s),",; \t");
		if(isDouble(s))
			*(double*)(prm)=atof(s);
		else fg=1;
		break;
	}
	case PRM_STRING:   	*(char**) (prm)=strdup(s);		break;
	case PRM_MAPIV:   	((MapIv *) (prm))->read(s);		break;
	case PRM_ENUM:		{
		int vl=readEnum(s);
		if(vl!=PRM_UNKNOWN) *(int *)(prm)=vl;
		else	fg=1;
		break;
	}
	default: fg=1; break;
	}
	return fg;
}
void Param::setVal(){
	switch(type){
	case PRM_INT:  		*(int *)  (prm)=value; break;
	case PRM_FG:		*(bool*)  (prm)=value; break;
	default: break;
	}
}

//===========================================================================================================
Param * readParam(char *name, char *value){
	Param *prm=findParam(name);
	if(prm==0) return 0;
	prm->readVal(value);
	return prm;
}

//===========================================================================================================
void Param::printDescr(){
	if(description==0) return;
	if(name ==0 || strlen(name)==0) {printf("\n====================== %s ====================== \n",description); return;}
	printf("-%s ", name);
	if(value==PRM_UNKNOWN) {
		if(type==PRM_INT) 	 printf("<int>");
		if(type==PRM_DOUBLE) printf("<float>");
		if(type==PRM_STRING) printf("<string>");
		if(type==PRM_FG)     printf("<0|1>");
	}
	if(enums){
		for(int i=0; enums[i]!=0; i++){
			if(i==0)printf("<");
			else    printf("|");
			printf("%s",enums[i]->name);
		}
		printf("> ");
	}
	printf("\t%s\n",description);
}
//============================================ Check if param name exists =============================
Param *findParam(const char * name){
	for(int i=0; pparams[i]!=0; i++){
		if(pparams[i]->name == 0) continue;
		if(keyCmp(pparams[i]->name,name)==0) return pparams[i];
	}
	return 0;
}

//============================================ Print Help page =========================================
void printHelp(){
	printf("\n");
	printf("The StereoGene program compares pairs of tracks and calculates kernel correlations\n");
	printf("Usage:\n");
	printf("$ ./StereoGene [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n");
	printf("\n");
	for(int i=0; pparams[i]!=0; i++){
		pparams[i]->printDescr();
		if(i%15==0 && i>0) {
			printf("Press q or Esc to exit or any key to go on");
			fflush(stdout);
			int c=xpause();
			if(c=='q' || c=='Q' || c==27) {printf("\n");break;}
			printf("\r");
		}
	}
}
//============================================ Read Config =========================================
void readConfig(char * cfg){
	FILE *f=gopen(cfg,"rt");

	if(f==0) return;
	char b[1024], *s;
	for(;(fgets(b,sizeof(b),f))!=0;){
		strtok(b,"\r\n");
		s=skipSpace(b);
		if(*s=='#') continue;
		if(*s==0) continue;
		s=skipSpace(b);
		if(*s==0) continue;
		readPrm(s);
	}
	fclose(f);

}
//============================================ Read Param =========================================
void readPrm(char *b){
	strtok(b,"#");
	char *s1=trim(strtok(b,"="));
	char *s2=trim(strtok(0,"="));
	readPrm(s1,s2);
}

void readPrm(char *key, char *val){
	Param* prm=findParam(key);
	if(prm!=0){
		if(prm->readVal(val)) errorExit("unknown value %s=%s",key,val);
	}
	else errorExit("unknown parameter %s",key);
}
//============================================ Parse comand line =========================================
void parseArgs(int argc, char **argv){
	char b[1024];
	//========================= Search for cgf ========================
	for(int i=1; i<argc; i++){
		if(*argv[i]=='-'){
			strcpy(b,argv[i]+1);
			if(keyCmp("cfg",b)==0) readConfig(argv[++i]);
		}
		if(strchr(argv[i],'=')){
			strcpy(b,argv[i]);
			char *s=strtok(b,"=");
			if(keyCmp("cfg",s)==0){
				readConfig(strtok(0,"= "));
			}
		}
	}
	//============================== Read command params =====================
	for(int i=1; i<argc; i++){
		if(*argv[i]=='-'){
			if(keyCmp("h",argv[i]+1)==0){
				printHelp();
				exit(0);
			}
			Param* prm=findParam(argv[i]+1);
			if(prm==0) errorExit("unknown parameter %s", argv[i]+1);
			if(prm->value!=PRM_UNKNOWN) prm->setVal();
			else if(prm->readVal(argv[++i])) errorExit("unknown value %s=%s",argv[i-1],argv[i]);
		}
		else if(strchr(argv[i],'=')){
			readPrm(argv[i]);
		}
		else{
			addFile(argv[i]);
		}
	}
	if(wStep==0)   wStep=wSize;
	if(RScriptFg) {writeDistCorr|=TOTAL; writeDistr=1;}
	if(complFg==0){
		if(strandFg0) complFg=COLLINEAR;
		else		  complFg=IGNORE_STRAND;
	}
	if(threshold < 1) threshold=1;
	profPath =makePath(profPath);
	trackPath=makePath(trackPath);
	resPath	 =makePath(resPath);
	if(verbose) silent=false;
}
//==============================================================================
//==============================================================================
//=========================================================== read value ===============================
// for debugging :
// set debugFg=DEBUG_LOG|DEBUG_PRINT
// set debS string for module identification//
// Use:  deb(n);  // print debug information as number
//       deb(format,...)    // print debug information as printf
//       deb(n,format,...)  // print debug information as a number and printf
// example:
// debS="fun1";
// deb(1);
// ....
// deb(2,"%i %f", n, d);
// ....
// deb("OK");

int main(int argc, char **argv) {
//	clearDeb();
//	debugFg=DEBUG_LOG|DEBUG_PRINT;
	const char * progName="StereoGene";
	verb("===== %s version %s =====\n",progName,version);
	char *chrom=getenv("SG_CHROM");
	if(chrom!=0) chromFile=strdup(chrom);

	unsigned long t=time(0);	id=t&0xffffff;	// define run id
	parseArgs(argc, argv);

	if(nfiles==0){
		printf("\n");
		printf("The %s program compares pairs of tracks and calculates kernel correlations\n",progName);
		printf("===========  version %s ========\n",version);
		printf("Usage:\n");
		printf("$ ./%s [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n",progName);
		printf("\n");
		printf("Say Stereogene -h for more information\n");
		printf("\n");
		exit(0);
	}
	if(aliaseFil!=0)  alTable.readTable(aliaseFil);		// read aliases
	readChromSizes(chromFile);							// read chromosomes

	if(cage) {CageMin(files[0].fname,files[1].fname); exit(0);}

	if(pcaFg) pcaMain(profile1);

	Correlator();

	if(logFile) {fclose(logFile); logFile=0;}	// debug log file
	return 0;
}

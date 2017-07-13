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
const int PRM_PATH=7;

const int PRM_UNKNOWN=-0XFFFFFFF;

//bool doLC=false;
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
	int value;					// a value that should be set if -prm is used in a command line
	int printFg;				// the param should be printed to PRM file
	const char *description;
	Param(const char* descr);
	Param(const char* _name, int print, int    *prm, const char* descr);
	Param(const char* _name, int print, int    *prm, int val, const char* descr);
	Param(const char* _name, int print, int    *prm, Name_Value **fg, const char* descr);
	Param(const char* _name, int print, bool   *prm, const char* descr);
	Param(const char* _name, int print, bool   *prm, bool val, const char* descr);
	Param(const char* _name, int print, double *prm, const char* descr);
	Param(const char* _name, int print, char * *prm, const char* descr);
	Param(const char* _name, int print, MapIv  *prm, const char* descr);
	Param(const char* _name, int print, char*  *_prm, const char* descr, bool path);

	void setVal();
	void init(const char* _name,int print, void* _prm, int type, Name_Value **fg, const char* descr);
	void printDescr();
	int readVal(char *s);
	int readEnum(char *s);
	char* printParamValue(char *buf);
};

const char* getNamebyVal(Name_Value **nval, int val){
	for(; nval!=0; nval++){
		if((*nval)->value == val) return (*nval)->name;
	}
	return "NA";
}

//===================================================================================================
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
Name_Value* LCScaleTypes[]={
		new Name_Value("LOG",LOG_SCALE),
		new Name_Value("LOG_LOG",LOG_LOG_SCALE),
		new Name_Value("LIN" ,LIN_SCALE),
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
		new Param("v"		    ,0, &verbose		,1, "verbose"),
		new Param("syntax"		,0, &syntax		,1, "strong syntax control in input files"),
		new Param("verbose"	    ,0, &verbose		,"verbose"),
		new Param("s"		    ,0, &silent		,1, "no output to stdout"),
		new Param("silent"	    ,0, &silent		,"no output to stdout"),


//=============================================================================================================
		new Param("preparation parameters"),
		new Param("bin" 	 	,1, &binSize  	,"bin size for input averaging"),
		new Param("scale" 	 	,0, &logScale  	,scaleTypes,"use Logarithmic or linear scale"),
		new Param("scaleFactor"	,0, &scaleFactor0 ,0),
		new Param("clear" 	  	,0, &clearProfile ,"force binary profile preparation"),
		new Param("c" 	  	    ,0, &clearProfile , 1,"force  binary profile preparation"),

//=============================================================================================================
		new Param("paths and files"),
		new Param("cfg" 		,0, &cfgFile 		,"config file"),
		new Param("profPath" 	,0, &profPath 	,"path for binary profiles", true),
		new Param("trackPath" 	,0, &trackPath 	,"path for tracks", true),
		new Param("resPath" 	,0, &resPath 		,"path for results", true),

		new Param("statistics" 	,0, &statFileName	 ,"cumulative file with statistics"),
		new Param("params" 		,0, &paramsFileName,"cumulative file with parameters"),
		new Param("log" 		,0, &logFileName	 ,"cumulative log-file"),
		new Param("aliases" 	,0, &aliaseFil	,0),
//=============================================================================================================
		new Param("input parameters"),
		new Param("chrom"		,0, &chromFile	,"chromosome file"),
		new Param("intervals"	,1, &intervFlag0	,intervalTypes,"interval type in BED file"),
		new Param("gene"		,0, &intervFlag0	,GENE		,"consider entire gene"),
		new Param("exon"		,0, &intervFlag0	,EXON		,"consider exons"),
		new Param("ivs"			,0, &intervFlag0	,IVS 	 	,"consider introns"),
		new Param("gene_beg"	,0, &intervFlag0	,GENE_BEG	,"gene starts"),
		new Param("exon_beg"	,0, &intervFlag0	,EXON_BEG	,"exons starts"),
		new Param("ivs_beg"		,0, &intervFlag0	,IVS_BEG	,"introns starts"),
		new Param("gene_end"	,0, &intervFlag0	,GENE_END	,"gene ends"),
		new Param("exon_end"	,0, &intervFlag0	,EXON_END	,"exons ends"),
		new Param("ivs_end"		,0, &intervFlag0	,IVS_END	,"introns ends"),
		new Param("bpType" 	  	,1, &bpType  		,bpTypes	,"The value used as a score for BroadPeak input file"),
		new Param("pcorProfile" ,1, &pcorProfile	,"Track for partial correlation"),
		new Param("NA"       	,1, &NAFlag     	,1 , "use NA values as unknown and fill them by noise"),
		new Param("threshold"	,1, &threshold	,"threshold for input data for removing too small values: 0..250"),
		new Param("map" 		,1, &mapFil 		,0),
		new Param("mapIv" 		,1, &miv			,0),
//=============================================================================================================
		new Param("Analysis parameters"),
		new Param("kernelType"	,1, &kernelType	,kernelTypes,0),
		new Param("kernelSigma"	,1, &kernelSigma  ,"Kernel width"),
		new Param("kernelShift"	,1, &kernelShift  ,0),
		new Param("wSize" 	  	,1, &wSize  		,"Window size"),
		new Param("wStep" 	  	,1, &wStep  		,0),
		new Param("kernelNS"	,0, &kernelNS  	,0),
		new Param("flankSize"	,1, &flankSize  	,0),
		new Param("maxNA"		,1, &maxNA  		,"Max number of NA values in window (percent)"),
		new Param("maxZero"		,1, &maxZero0  	,"Max number of zero values in window (percent)"),
		new Param("nShuffle"	,1, &nShuffle  	,"Number of shuffles for background calculation (percent of window pairs)"),
		new Param("MaxShuffle"	,0, &maxShuffle  	,"Max number of shuffles"),
		new Param("MinShuffle"	,0, &minShuffle  	,"Min number of shuffles"),
		new Param("noiseLevel"	,1, &noiseLevel  	,0),
		new Param("complFg"		,1, &complFg		,complFlags,0),
//=============================================================================================================
		new Param("Output parameters"),
		new Param("outSpectr" 	,0, &outSpectr   	,"write fourier spectrums"),
		new Param("outChrom" 	,1, &outChrom   	,"write statistics by chromosomes"),
		new Param("writeDistr" 	,1, &writeDistr   ,"write foreground and background distributions"),
		new Param("Rscrpit" 	,0, &RScriptFg   	,0),
		new Param("r" 			,0, &RScriptFg   	,1,"write R script for the result presentation"),
		new Param("crossWidth" 	,0, &crossWidth   ,0,"Width of cross-correlation plot"),
		new Param("Distances" 	,1, &writeDistCorr,1,"Write distance correlations"),
		new Param("outLC"		,1, &outWIG		,outWigTypes,"parameters for local correlation file"),
		new Param("lc"			,0, &outWIG		,WIG_BASE|WIG_SUM,"produce profile correlation with parameter BASE"),
		new Param("LCScale"		,1, &LCScale		,LCScaleTypes,"Local correlation scale: LOG_LOG | LOG | LIN"),

		new Param("outThreshold",1, &outThreshold	,"threshold for output to correlation profile scaled to 0..1000"),
		new Param("corrOnly" 	,0, &corrOnly   	,0),
		new Param("corr" 		,0, &corrOnly   	,1, 0),
		new Param("outRes" 		,0, &outRes 		,outResTypes,"format for results in statistics file"),
		new Param("AutoCorr"  	,0, &doAutoCorr  	,0),
//=========================================== Additional parameters (see Undocumented) ===============================
		new Param("outBPeak" 	,0, &writeBPeak   ,"write BradPeak file"),
		new Param("pVal"		,0, &pVal  		,"threshold for BroadPeak output"),
		new Param("qVal"		,0, &qVal  		,"threshold for BroadPeak output"),
		new Param("pcaSegment" 	,0, &pcaSegment	,0),
		new Param("nPca" 		,0, &nPca			,0),
		new Param("cage" 		,0, &cage			,0),
		new Param("inpThreshold",0, &inpThreshold ,0),	//input binarization testing, %of max
		new Param("d"			,0, &debugFg  	,1, 0),	//debug mode
		new Param("pdf"			,0, &writePDF  ,1, 0),	//write R plots to pdf

		new Param("Happy correlations!"),
		0,
};

//=================================================================================================
//===================================  End declaration ============================================
//=================================================================================================
void Param::init(const char* _name,int print,void *_prm, int _type, Name_Value **fg, const char* descr){
	name=_name;
	printFg=print;
	enums=fg;
	description=descr;
	prm=_prm;
	type=_type;
	value=PRM_UNKNOWN;
}
Param::Param(const char* descr)
	{init("",0,0,0,0 ,descr);}
Param::Param(const char* _name,int print,  int    *_prm, const char* descr)
	{init(_name,print,_prm,PRM_INT	,0 ,descr);}
Param::Param(const char* _name,int print,  int    *_prm, Name_Value **fg, const char* descr)
	{init(_name,print,_prm,PRM_ENUM	,fg,descr);}
Param::Param(const char* _name,int print,  bool   *_prm,  const char* descr)
	{init(_name,print,_prm,PRM_FG		,0 ,descr);}
Param::Param(const char* _name,int print,  double *_prm,  const char* descr)
	{init(_name,print,_prm,PRM_DOUBLE	,0 ,descr);}
Param::Param(const char* _name,int print,  char*  *_prm, const char* descr)
	{init(_name, print,_prm,PRM_STRING	,0 ,descr);}
Param::Param(const char* _name,int print,  char*  *_prm, const char* descr, bool path)
	{init(_name,print,_prm,PRM_PATH	,0 ,descr);}
Param::Param(const char* _name,int print, MapIv   *_prm, const char* descr)
	{init(_name,print,_prm,PRM_MAPIV	,0 ,descr);}
Param::Param(const char* _name,int print,  int    *_prm, int val, const char* descr)
	{init(_name,print,_prm,PRM_INT	,0 ,descr);	value=val;}
Param::Param(const char* _name,int print,  bool    *_prm, bool val, const char* descr)
	{init(_name,print,_prm,PRM_FG	,0 ,descr);	value=val;}
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
	case PRM_STRING:
		if(s==0 || strlen(s)) *(char**) (prm)=strdup(s);
		else prm=0;
		break;
	case PRM_PATH:
		if(s==0 || strlen(s)) *(char**) (prm)=makePath(s);
		else prm=0;
		break;
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
char *Param::printParamValue(char *buf){
	strcpy(buf,"NONE");
	switch(type){
	case PRM_INT: 		sprintf(buf,"%i",*(int *)prm); break;
	case PRM_DOUBLE: 	sprintf(buf,"%f",*(double*)prm); break;
	case PRM_STRING: 	if(prm){
		char *s=*(char**)prm;
		if(s) sprintf(buf,"%s",s);} break;
	case PRM_ENUM: 		sprintf(buf,"%s",getNamebyVal(enums,*(int*)prm)); break;
	case PRM_FG: 		sprintf(buf,"%i",(*(int*)prm) ? 1:0); break;
	case PRM_MAPIV:		if(prm)((MapIv *)prm) -> print(buf); break;
	case PRM_PATH:		if(prm){
		char *s=*(char**)prm;
		if(s) sprintf(buf,"%s",s);} break;
	}
	return buf;
}

void printParamNames(FILE* f){
	fprintf(f,"id");
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"\t%s",pparams[i]->name);
	}
}
void printParams(FILE* f){
	char b[256];
	fprintf(f,"%08lx",id);
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"\t%s",pparams[i]->printParamValue(b));
	}
}
void printXMLparams(FILE *f){
	char b[256];
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"%s=\"%s\" ",pparams[i]->name,pparams[i]->printParamValue(b));
	}
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
		if(type==PRM_PATH) printf("<string>");
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
	if(keyCmp(key,"in")==0) {addFile(val); return;}
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
	if(RScriptFg) {writeDistCorr=1; writeDistr=1;}
	if(complFg==0){complFg=IGNORE_STRAND;}
	if(threshold < 1) threshold=1;
	profPath =makePath(profPath);
	trackPath=makePath(trackPath);
	resPath	 =makePath(resPath);
	if(verbose) silent=false;
}
//==============================================================================
//==============================================================================
//================================================== search appropriate cfg file
//void readCfg(int argc, const char *argv[]) {
//	argv[0]=correctFname(argv[0]);
//	char b[1024];
//  getFnameWithoutExt(b,argv[0]);
//	strcat(b,".cfg");
//	char *cfg=cfgName(b, (char*)"cfg");
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
//==============================================================================
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
//	test();
//	clearDeb();
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	for(int i=0; i<argc; i++){strtok(argv[i],"\r\n");}

	const char * progName="StereoGene";
	char *chrom=getenv("SG_CHROM");
	if(chrom!=0) chromFile=strdup(chrom);
	unsigned long t=time(0);	id=(unsigned int)t;	// define run id
	parseArgs(argc, argv);
	if(debugFg){
		clearDeb();
		debugFg=DEBUG_LOG|DEBUG_PRINT;
	}
	makeDirs();
	if(strcmp(logFileName,"null")==0 || strcmp(logFileName,"NULL")==0) logFileName=0;
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
//	if(cage) {CageMin(files[0].fname,files[1].fname); exit(0);}
//PrintParams();
//exit(0);
	Correlator();
	fflush(stdout);
//	fclose(stdout);
	return 0;
}

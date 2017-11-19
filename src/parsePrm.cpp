/*
 * parsePrm.cpp
 *
 *  Created on: Feb 2, 2017
 *      Author: mironov
 */
#include "track_util.h"
#include <unistd.h>

#if defined(_WIN32)
#include <conio.h>
int nHelpLines=20;
#else
#include <termios.h>
int nHelpLines=2000;
#endif

int xpause(){
	int c=0;
#if defined(_WIN32)
	c=_getch(); printf("\n%c\n",c);
#else
//	return 0;
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
const int PRM_PATH=7;

const int PRM_UNKNOWN=-0XFFFFFFF;


struct Param{
	const char* name;			// command line (cfg) argument name
	int type;					// parameter type
	Name_Value **enums;			// array of aviable values
	void *prm;					// pointer to the argument value
	int value;					// a value that should be set if -prm is used in a command line
	int printFg;				// the param should be printed to PRM file
	int prog;					// the parameter relevant to the program of given type
	const char *description;
	Param(int prg, const char* descr);
	Param(int prg, const char* _name, int print, int    *prm, const char* descr);
	Param(int prg, const char* _name, int print, int    *prm, int val, const char* descr);
	Param(int prg, const char* _name, int print, int    *prm, Name_Value **fg, const char* descr);
	Param(int prg, const char* _name, int print, bool   *prm, const char* descr);
	Param(int prg, const char* _name, int print, bool   *prm, bool val, const char* descr);
	Param(int prg, const char* _name, int print, double *prm, const char* descr);
	Param(int prg, const char* _name, int print, char * *prm, const char* descr);
	Param(int prg, const char* _name, int print, char*  *_prm, const char* descr, bool path);

	void setVal();
	void init(int prg, const char* _name,int print, void* _prm, int type, Name_Value **fg, const char* descr);
	void printDescr();
	int readVal(char *s);
	int readEnum(char *s);
	char* printParamValue(char *buf);
	char* printParamXML(char *buf);
};

const char* getNamebyVal(Name_Value **nval, int val){
	for(; nval!=0; nval++){
		if((*nval)->value == val) return (*nval)->name;
	}
	return "NA";
}

//===================================================================================================
Param *findParam(const char * name);
void readPrm(char *s);
void readPrm(char *key, char *val);
void printHelp();


//==============================================================================
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
Name_Value* LCFlags[]={
		new Name_Value("BASE",BASE),
		new Name_Value("CENTER",CENTER),
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
		new Param( 7,"common parameters"),
		new Param( 7,"v"		    ,0, &verbose	,1, "verbose"),
		new Param( 7,"syntax"		,0, &syntax		,1, "strong syntax control in input files"),
		new Param( 7,"verbose"		,0, &verbose	,   "verbose"),
		new Param( 7,"s"		    ,0, &silent		,1, "no output to stdout"),
		new Param( 7,"silent"		,0, &silent		,	"no output to stdout"),
//======================== =====================================================================================
		new Param( 7,"preparation parameters"),
		new Param( 7,"bin" 	 	,1, &binSize  		,"bin size for input averaging"),
		new Param( 7,"clear" 	,0, &clearProfile 	,"force binary profile preparation"),
		new Param( 7,"c" 	  	,0, &clearProfile 	, 1,"force  binary profile preparation"),
//======================== =====================================================================================
		new Param( 7,"paths and files"),
		new Param( 7,"cfg" 			,0, &cfgFile 		,"config file"),
		new Param( 7,"profPath" 	,1, &profPath 		,"path for binary profiles", true),
		new Param( 7,"trackPath" 	,1, &trackPath 		,"path for tracks", true),
		new Param(SG,"resPath" 		,1, &resPath 		,"path for results", true),
		new Param( 7,"confounder"	,0, &confFile 		,"confounder filename"),

		new Param(SG,"statistics"	,0, &statFileName	,"cumulative file with statistics"),
		new Param(SG,"params" 		,0, &paramsFileName	,"cumulative file with parameters"),
		new Param( 7,"log" 		,0, &logFileName	,"cumulative log-file"),
		new Param( 7,"aliases" 	,0, &aliaseFil		,0),
		new Param( 7,"id_suff" 	,0, &idSuff		,0),
//======================== =====================================================================================
		new Param( 7, "input parameters"),
		new Param( 7, "chrom"		,0, &chromFile	,"chromosome file"),
		new Param( 7, "BufSize"	,0, &binBufSize	,"Buffer Size"),
		new Param( 7, "bpType" 		,1, &bpType  		,bpTypes	,"The value used as a score for BroadPeak input file"),
		new Param(SG|PRJ,"pcorProfile" ,1, &pcorProfile	,"Track for partial correlation"),
		new Param(PRJ,"outPrjBGr"   ,0, &outPrjBGr		,"Write BedGraph for projections"),
		new Param(SG, "NA"       	,1, &NAFlag     	,1 , "use NA values as unknown and fill them by noise"),
		new Param(SG, "threshold"	,1, &threshold	,"threshold for input data for removing too small values: 0..250"),
//======================== =====================================================================================
		new Param(SG, "Analysis parameters"),
		new Param(SG, "kernelType"	,1, &kernelType	,kernelTypes,0),
		new Param(SG, "customKern"	,1, &customKern	,0),
		new Param(SG, "kernelSigma"	,1, &kernelSigma,"Kernel width"),
		new Param(SG, "kernelShift"	,1, &kernelShift,0),
		new Param(SG, "wSize" 	  	,1, &wSize  	,"Window size"),
		new Param(SG, "wStep" 	  	,1, &wStep  	,0),
		new Param(SG, "kernelNS"	,0, &kernelNS  	,0),
		new Param(SG, "flankSize"	,1, &flankSize  ,0),
		new Param(SG, "maxNA"		,1, &maxNA0  	,"Max number of NA values in window (percent)"),
		new Param(SG, "maxZero"		,1, &maxZero0  	,"Max number of zero values in window (percent)"),
		new Param(SG, "nShuffle"	,1, &nShuffle  	,"Number of shuffles for background calculation"),
		new Param(SG, "noiseLevel"	,1, &noiseLevel ,0),
		new Param(SG, "complFg"		,1, &complFg	,complFlags,0),
		new Param(SG, "LCFg"		,1, &lcFlag		,LCFlags,0),
//======================== =====================================================================================
		new Param(SG, "Output parameters"),
		new Param(SG, "outSpectr" 	,1, &outSpectr    ,"write fourier spectrums"),
		new Param(SG, "outChrom" 	,1, &outChrom     ,"write statistics by chromosomes"),
		new Param(SG, "writeDistr" 	,1, &writeDistr   ,"write foreground and background distributions"),
		new Param(SG, "Rscript" 	,1, &RScriptFg    ,0),
		new Param(SG, "r" 			,0, &RScriptFg    ,1,"write R script for the result presentation"),
		new Param(SG, "crossWidth" 	,0, &crossWidth   ,0,"Width of cross-correlation plot"),
		new Param(SG, "Distances" 	,1, &writeDistCorr,1,"Write distance correlations"),
		new Param(SG, "outLC"		,1, &outLC		  ,"parameters for local correlation file"),
		new Param(SG, "lc"			,0, &outLC		  ,1,"produce profile correlation"),
		new Param(SG, "LCScale"		,0, &LCScale	  ,LCScaleTypes,"Local correlation scale: LOG | LIN"),
		new Param(SG, "L_FDR"		,1, &LlcFDR	      ,"threshold on left FDR when write the local correlation"),
		new Param(SG, "R_FDR"		,1, &RlcFDR	      ,"threshold on right FDR when write the local correlation"),
		new Param(SG, "outRes" 		,0, &outRes 	  ,outResTypes,"format for results in statistics file"),
		new Param(SG, "AutoCorr"  	,1, &doAutoCorr   ,0),
//======================== =================== Additional parameters (see Undocumented) ===============================
		new Param(SG, "inpThreshold",0, &inpThreshold ,0),	//input binarization testing, %of max
		new Param(7 , "debug"		,0, &debugFg   	  ,0),	//debug mode
		new Param(7 , "d"			,0, &debugFg   	  ,1, 0),	//debug mode
		new Param(SG, "pdf"			,1, &writePDF  	  ,1, 0),	//write R plots to pdf

		new Param(PG, "pgLevel"		,1, &pgLevel  	  ,1, 0),	//minimal level in ENCODE to be taken into account

		new Param(7, "Happy correlations!"),
		0,
};

//=================================================================================================
//===================================  End declaration ============================================
//=================================================================================================
void Param::init(int prg, const char* _name,int print,void *_prm, int _type, Name_Value **fg, const char* descr){
	name=_name;
	printFg=print;
	enums=fg;
	description=descr;
	prm=_prm;
	type=_type;
	value=PRM_UNKNOWN;
	prog=prg;
}
Param::Param(int prg, const char* descr)
	{init(prg,"",0,0,0,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  int    *_prm, const char* descr)
	{init(prg,_name,print,_prm,PRM_INT	,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  int    *_prm, Name_Value **fg, const char* descr)
	{init(prg,_name,print,_prm,PRM_ENUM	,fg,descr);}
Param::Param(int prg, const char* _name,int print,  bool   *_prm,  const char* descr)
	{init(prg,_name,print,_prm,PRM_FG		,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  double *_prm,  const char* descr)
	{init(prg,_name,print,_prm,PRM_DOUBLE	,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  char*  *_prm, const char* descr)
	{init(prg,_name, print,_prm,PRM_STRING	,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  char*  *_prm, const char* descr, bool path)
	{init(prg,_name,print,_prm,PRM_PATH	,0 ,descr);}
Param::Param(int prg, const char* _name,int print,  int    *_prm, int val, const char* descr)
	{init(prg,_name,print,_prm,PRM_INT	,0 ,descr);	value=val;}
Param::Param(int prg, const char* _name,int print,  bool    *_prm, bool val, const char* descr)
	{init(prg,_name,print,_prm,PRM_FG	,0 ,descr);	value=val;}
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
	case PRM_DOUBLE: 	sprintf(buf,"%.2g",*(double*)prm); break;
	case PRM_STRING: 	if(prm){
		char *s=*(char**)prm;
		if(s) sprintf(buf,"%s",s);} break;
	case PRM_ENUM: 		sprintf(buf,"%s",getNamebyVal(enums,*(int*)prm)); break;
	case PRM_FG: 		sprintf(buf,"%i",(*(int*)prm) ? 1:0); break;
	case PRM_PATH:		if(prm){
		char *s=*(char**)prm;
		if(s) sprintf(buf,"%s",s);} break;
	}
	return buf;
}

char *Param::printParamXML(char *buf){
	char bb[1024];
	sprintf(buf,"%s=\"%s\"",name,printParamValue(bb));
	return buf;
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
	char *prm=strtok(b,"#");
	char *val=strchr(b,'=');
	if(*val==0)  val=b+strlen(b);
	else 		*val++=0;
	readPrm(trim(prm), trim(val));
}

void readPrm(char *key, char *val){
	if(keyCmp(key,"in")==0) {addFile(val); return;}
	if(keyCmp(key,"cfg")==0) readConfig(val);
	Param* prm=findParam(key);
	if(prm!=0){
		if(prm->readVal(val)) errorExit("unknown value %s=%s",key,val);
	}
	else errorExit("unknown parameter {%s}",key);
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
//	if(threshold < 1) threshold=1;
	if(customKern) kernelType=KERN_CUSTOM;
	profPath =makePath(profPath);
	trackPath=makePath(trackPath);
	resPath	 =makePath(resPath);
	if(verbose) silent=false;
}

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
//==============================================================================
void printParamNames(FILE* f){
	fprintf(f,"id\tversion");
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"\t%s",pparams[i]->name);
	}
	fprintf(f,"\n");
}
void printParams(FILE* f){
	char b[256];
	fprintf(f,"%08lx\t%s",id,version);
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"\t%s",pparams[i]->printParamValue(b));
	}
	fprintf(f,"\n");
}
void printXMLparams(FILE *f){
	char b[256];
	for(int i=0; pparams[i] ; i++){
		if(pparams[i]->printFg) fprintf(f,"%s=\"%s\" ",pparams[i]->name,pparams[i]->printParamValue(b));
	}
}


void printHelp(){
	printProgDescr();
	for(int i=0, j=0; pparams[i]!=0; i++){
		if(!(pparams[i]->prog & progType)) continue;
		pparams[i]->printDescr();
		if(j%nHelpLines==0 && j>0) {
			printf("Press q or Esc to exit or any key to go on");
			fflush(stdout);
			int c=xpause();
			printf("\n");
			if(c=='q' || c=='Q' || c==27) {break;}
		}
		 j++;
	}
}


void initSG(int argc, char **argv){
	for(int i=0; i<argc; i++){strtok(argv[i],"\r\n");}

	char *chrom=getenv("SG_CHROM");
	if(chrom!=0) chromFile=strdup(chrom);
	unsigned long t=time(0);	id=(unsigned int)t;	// define run id
	parseArgs(argc, argv);
	makeDirs();
	if(strcmp(logFileName,"null")==0 || strcmp(logFileName,"NULL")==0) logFileName=0;
	if(nfiles==0) printMiniHelp();
	if(aliaseFil!=0)  alTable.readTable(aliaseFil);		// read aliases
	readChromSizes(chromFile);							// read chromosomes
}






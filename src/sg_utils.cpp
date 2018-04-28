/*
 * sg_utils.cpp
 *
 *  Created on: 05 дек. 2017 г.
 *      Author: andrey
 */

#include "track_util.h"
unsigned int hashx(unsigned int h,char c){
	return h+(c-32)+1234;
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
	id=0;
	id=hashx(id,chromFile);
	id=hashx(id,flankSize);
	id=hashx(id,kernelType);
	id=hashx(id,binSize);
	id=hashx(id,nShuffle);
	id=hashx(id,maxZero);
	id=hashx(id,maxNA0);
	id=hashx(id,kernelSigma);
	id=hashx(id,wSize);
	id=hashx(id,wStep);
	id=hashx(id,threshold);
	id=hashx(id,outFile);
	id=hashx(id,mtime());
}

//====================================================================================
ScoredRange::ScoredRange(){
	chr=0; chrom=0;	beg=end=0;	score=0;
}
//==================================== print a range to a BED GRAPH
void ScoredRange::printBGraph(FILE *f){
	if(score!=NA){
		fprintf(f,"%s\t%li\t%li\t",chrom,beg, end);
		double xx=abs(score);
		if(xx < 0.000001) 	fprintf(f,"0.0\n");
		else if(xx <= 0.0000999) 	fprintf(f,"%.7f\n",score);
		else if(xx <= 0.000999) fprintf(f,"%.6f\n",score);
		else if(xx <= 0.00999)  fprintf(f,"%.5f\n",score);
		else if(xx <= 0.0999)   fprintf(f,"%.5f\n",score);
		else                fprintf(f,"%.3f\n",score);
	}
}


//====================================================================================
//===============================================  Chromosomes
Chromosome::Chromosome(char *chr, long l, int bb){
    chrom=chr;		//Chromosome name
    length=l;		//Chromosome length
    base=bb;		//Start position in the binary profile
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


int CHROM_BUFF=300;

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
		s=trim(buff);
		if(*s==0 || *s=='#') 	continue;

		if(n_chrom>=max_chrom) {
			max_chrom+=CHROM_BUFF;
			chrom_list=(Chromosome *)realloc(chrom_list,max_chrom*sizeof(Chromosome));
		}
	    if((s1=strtok(s," \t\r\n"))==NULL) continue;
	    if((s2=strtok(0," \t\r\n"))==NULL) continue;

	    chrom_list[n_chrom]=Chromosome(strdup(s1), atol(s2), profileLength);

		int filLen=(chrom_list[n_chrom].length+binSize-1)/binSize;
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
	long p=pos/binSize+curChrom->base;
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

//========================================================================================
void filePos2Pos(int pos, ScoredRange *gr, int length){
	Chromosome *ch0=getChromByPos(pos);

	if(ch0==0) return;
	pos-=ch0->base;
	long p1=binSize*long(pos);
	gr->chr=ch0;
	gr->chrom=ch0->chrom;
	gr->end=(gr->beg=p1)+length;
	if(gr->end >= ch0->length) gr->end = ch0->length-1;
	return;
}

//========================================================================================
int inputErr;		// flag: if input track has errors
int inputErrLine;	// Error line in the input
char curFname[4048];	// current input file

//========================================================================================
Chromosome *checkRange(ScoredRange *gr){
	Chromosome* chr=findChrom(gr->chrom);
	if(chr==0) return 0;

	if(gr->beg < 0){
		if(inputErr == 0){ inputErr=1;
			writeLogErr(
					"File <%s> line #%d: incorrect segment start: chrom=%s  beg=%ld.  Ignored\n",curFname, inputErrLine,gr->chrom, gr->beg);
		}
		return 0;
	}
	if(gr->end  >= chr->length){
		if(inputErr == 0){ inputErr=1;
			writeLogErr(
					"File <%s> line #%d: incorrect segment end: chrom=%s  end=%ld.  Ignored\n",curFname, inputErrLine,gr->chrom, gr->end);
		}
		return 0;
	}
	if(gr->end < gr->beg) return 0;
	return chr;
}

//======================================================================
BufFile::~BufFile(){
	if(f!=0) fclose(f);
	if(buffer) xfree(buffer,"buff file");
}

void BufFile::init(const char *fname){
	f=fopen(fname,"rb"); buffer=0;
	getMem0(buffer,SG_BUFSIZ+SG_BUFEXT,"err");
	int n=fread(buffer,1,SG_BUFSIZ,f);
	if(n <= 0) {curString=0;}
	else {curString=buffer; buffer[n]=0;}
}


char *BufFile::getString(){
	if(curString==0) return 0;
	char *s0=curString;
	while(isspace(*s0)) s0++;
	char *ss=strchr(curString,'\n');
	if(ss==0){
		strcpy(buffer,curString);
		char *bb=buffer+strlen(curString);
		int n=fread(bb,1,SG_BUFSIZ,f);
		if(n <= 0) {curString=0; return s0;}
		else {curString=buffer; bb[n]=0;}
		return(getString());
	}
	curString=ss+1;
	while(isspace(*ss)) *ss--=0;
//	*ss=0;
	return s0;
}

// get major version number (1.64.2 -> 1.64)
char * getMajorVer(const char *ver, char *buf){
	char b[1024];
	strcpy(b,ver);
	strcpy(buf,strtok(b,"."));
	strcat(buf,".");
	strcat(buf,strtok(0,"."));
	return buf;
}
//========== Convert kernel type to a string
char kernType[1024];
const char*getKernelType(){
	const char *type;
	if(kernelType==KERN_NORM	 ) type="N";
	else if(kernelType==KERN_LEFT_EXP ) type="L";
	else if(kernelType==KERN_RIGHT_EXP) type="R";
	else if(kernelType==KERN_CUSTOM) type=customKern;
	else return "X";
	char b[80];
	if(kernelShift >0){
		sprintf(b,"%s_%.1fR",type,kernelShift/1000);
		return strcpy(kernType,b);
	}
	else if(kernelShift <0){
		sprintf(b,"%s_%.1fL",type,-kernelShift/1000);
		return strcpy(kernType,b);
	}
	else return type;
}


//================== make path - add '/' if necessary
char* makePath(char* pt){
	if(pt==0) return pt;
	char b[2048];
	char *s=pt+strlen(pt)-1;
	if(*s=='/') *s=0;
	return strdup(strcat(strcpy(b,pt),"/"));
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

//===================== convert text flag to a binary
int getFlag(char*s){
	int fg=0;
	if(		keyCmp(s,"1")==0 || keyCmp(s,"YES")==0 || keyCmp(s,"ON" )==0) {fg=1;}
	else if(keyCmp(s,"0")==0 || keyCmp(s,"NO")==0  || keyCmp(s,"OFF")==0) {fg=0;}
	else fg=-1;
	return fg;
}
//=================================================================
//============================= File list =========================
//=================================================================
int   fileId=0;

void addFile(char* fname, int id){
	fname=trim(fname);
	if(strlen(fname)==0) return;

	files[nfiles].fname=strdup(fname);
	files[nfiles].listId=fileId;
	nfiles++;
}

void addFile(char* fname){
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
			s=trim(b);
			if(strlen(s)==0 || *s=='#') continue;
			addFile(s, fileId);
		}
		fclose(f); fileId++;
		return;
	}
	else{
		addFile(fname, fileId); fileId++;
	}
}




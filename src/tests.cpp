/*
 * tests.cpp
 *
 *  Created on: 22.02.2013
 *      Author: 1
 */
#include <stdio.h>


#include "track_util.h"

int wigStep=300;
void printBytes(unsigned char* byteprofile,int from, int to){
	printf("==========\n");
	for(int i=from; i<to; i+=20){
		printf("%4i\t",i);
		for(int j=0; j<20 && i+j< to; j++)
			if(byteprofile[i+j]==0) printf("  NA");
			else printf(" %3i",byteprofile[i+j]);
		printf("\n");
	}
}


double *readDistr(const char *fname, int &nn){
	int capacity=1000;
	char b[256];
	double *dstr; getMem(dstr,capacity, "readDistr");
	nn=0;
	FILE *f=fopen(fname,"rt");
	for(;fgets(b,sizeof(b),f);){
		strtok(b," \t\n\r");
		double a=atof(b);
		if(nn>=capacity) capacity*=2;
		dstr=(double*)realloc(dstr,capacity*sizeof(double));
		dstr[nn++]=a;
	}
	fclose(f);
	return dstr;
}

double *gauss(int l){
	double *x; getMem(x,l, "gauss");
	double a=0.99;
	x[0]=0;
	for(int i=1; i<l; i++)
		x[i]=a*x[i-1]+rGauss(0,0.1);
	return x;
}


void factor(int n){
	int k=sqrt(n)+1;
	for(int i=2; i<k; i++){
		while(n%i ==0){
			printf(" %i;",i);
			n/=i;
		}
	}
	if(n!=1) printf(" %i;",n);
}

//============================================================== print profile (testing)
void printProfile(float *profile, int from, int to){
	for(int i=from; i< to; i+=20){
		printf("%i\t",i);
		for(int j=0; j<20 && i+j< to; j++){
			if(profile[i+j]==NA) printf("    NA");
			else printf(" %5.1f",profile[i+j]*100);
		}
		printf("\n");
	}
}
//============================================================== print profile (testing)
void printProfile(float *profile){
	printProfile(profile, 0, profileLength);
}
//=========================================================== print chromosome data (testing)
void printChrom(){
   for(int i=0; i<n_chrom; i++){
      Chromosome chrom=chrom_list[i];
      printf("%s\t%li\t%i\n",chrom.chrom,chrom.length,chrom.base);
   }
}

//===========================================================
void bTrack::printBytes(int from, int to){
	printBytes(stdout, from,to);
}
void bTrack::printBytes(FILE *f, int from, int to){
	fprintf(f,">>>>>>>>> %i..%i\n",from,to);
	for(int i=from; i<to; i++) {if(i%25==0) fprintf(f,"\n");fprintf(f," %3i",bytes[i]);}
	fprintf(f,"\n");
	if(hasCompl){
		fprintf(f,"<<<<<<<<<");
		for(int i=from; i<to; i++)  {if(i%25==0) fprintf(f,"\n");fprintf(f," %3i",cbytes[i]);}
		fprintf(f,"\n");
	}
}
void bTrack::writeBytes(FILE*f){ // Print bytes (for testing)
	int n=0, t=15;
	for(int i=0; i<lProf; i++) {
		if(bytes[i]>15) {fprintf(f," %3i\n",bytes[i]); n++;}
	}
	printf("t=%i n=%i\n",t,n);
}

void bTrack::printWindow(int n){
	printWindow(stdout,n);
}
void bTrack::printWindow(FILE *f,int n){
	fprintf(f,"==========\n");
	for(int i=0; i<n; i++) fprintf(f," %4.1f",profWindow[i]);
	fprintf(f,"\n");
}
//===============================================================================
//===============================================================================
//===============================================================================
struct GenePoint{
	char *chrom;
	int   pos;
	char  strand;
	GenePoint();
	GenePoint(const char*s, int p, char strnd);
};

GenePoint::GenePoint(){
	chrom=0; pos=0; strand='+';
}

GenePoint::GenePoint(const char*s, int p, char strnd){
	chrom=strdup(s);
	pos=p;
	strand=strnd;
}

GenePoint **gPoints;
int n_points=0;
int maxPoints=20000;
float *dens;
int tstWindow=10000;

int GPCmp(const void *o1, const void *o2){
	GenePoint **gp1=(GenePoint**) o1, **gp2=(GenePoint**) o2;
	GenePoint *g1=*gp1, *g2=*gp2;
	int a=strcmp(g1->chrom,g2->chrom);
	if(a) return a;
	return g1->pos-g2->pos;
}

int chkPos(const char *chr, int pos, int i){
	int fg=strcmp(chr,gPoints[i]->chrom);
	if(fg) return fg;
	if(pos < gPoints[i]->pos-tstWindow) return -1;
	if(pos > gPoints[i]->pos+tstWindow) return  1;
	return 0;
}

int findPoint(const char* chr, int pos){
	int i1=0, i2=n_points,i=0,fg=0;

	while(i2-i1 > 1){
		i=(i1+i2)/2;
		fg=chkPos(chr,pos,i);
		if(fg==0) return i;
		if(fg<0) i2=i;
		if(fg>0) i1=i;
	}
	if(chkPos(chr,pos,i)==0) return i;
	return -1;
}

void readPoints(const char *fname){
	char b[2096],*s, chrBuf[80];
	makeFileName(b,trackPath,fname);
	FILE *f=xopen(b,"rt");
	getMem(gPoints,maxPoints, "readPoints");
	char strand='+';
	for(; (s=fgets(b,sizeof(b),f))!=0;){
		strtok(s,"\n\r");									 	//======== remove end-of-line signes
		if(*s=='#') continue;
		s=strtok(b,"\t "); strcpy(chrBuf,s);
		int p=atoi(s=strtok(0,"\t "));					//======== read score
		s=strtok(0," \t\n"); if(s==0) break;			//======== skip the second position
		s=strtok(0," \t\n"); if(s==0) break;			//======== skip the name
		s=strtok(0,"\t ");								//======== skip the score
		s=strtok(0,"\t "); strand=*s;
		if(strand!='+') continue;

		gPoints[n_points++]=new GenePoint(chrBuf,p,strand);
		if(n_points >= maxPoints) break;
	}
	fclose(f);
	qsort(gPoints, n_points, sizeof(GenePoint *), GPCmp);
}

void addSgm(int gPos, int beg, int end, int score, char strand){
//if(beg > 925500) deb(">>#0  gPos=%i  beg=%i..%i xpos=%i score=%i",gPos,beg,end,beg-gPos,score);
	for(int i=beg; i<end; i++){
		int pp=i-gPos;
		if(strand=='-') pp=-pp;
		if(pp< -tstWindow || pp>=tstWindow) continue;
		int p=pp+tstWindow;
		dens[p]+=score;
	}
}

void addSegment(char *chr, int beg, int end, int score){
	int ii=findPoint(chr,beg);
	if(ii<0) return;
	for(int i=ii; i>=0; i--){
		if(strcmp(gPoints[i]->chrom,chr)) break;
		if(beg > gPoints[i]->pos + tstWindow) break;
		addSgm(gPoints[i]->pos,beg,end,score, gPoints[i]->strand);
	}
	for(int i=ii+1; i<n_points; i++){
		if(strcmp(gPoints[i]->chrom,chr)) break;
		if(end < gPoints[i]->pos - tstWindow) break;
		addSgm(gPoints[i]->pos,beg,end,score, gPoints[i]->strand);
	}
}

void readWig(const char *fname){
	char b[2096],abuf[1000], chBuf[100],  buff[4096],*s, *chrom=(char*)"";
	int fg=0, span=0, start=0, beg=0, end=0, score=0, step;
	makeFileName(b,trackPath,fname);
	FILE *f=xopen(b,"rt");
	int ii=0;

	for(; (s=fgets(buff,sizeof(buff),f))!=0; ii++){
		if(ii%1000000 == 0)
			verb("%ik\t%s\n",ii/1000,chrom);
		strtok(s,"\n\r");									 	//======== remove end-of-line signes
		if(*s=='#') continue;
		if(strncmp(s,"variableStep", 12)==0){			//======== read parameters for variableStep
			chrom=getAttr(buff,(char *)"chrom",chBuf);
			span=1; fg=0;
			s=getAttr(buff,(char *)"span",abuf);
			if(s!=0) span=atoi(s);
			continue;
		}
		else if(strncmp(s,(char *)"fixedStep", 9)==0){	//======== read parameters for fixedStep
			fg=1;
			chrom=getAttr(buff,(char *)"chrom",chBuf);
			s=getAttr(buff,(char *)"span",abuf);
			if(s!=0) span=atoi(s);
			s=getAttr(buff,(char *)"step",abuf);
			if(s!=0) step=atoi(s);
			s=getAttr(buff,(char *)"start",abuf);
			if(s!=0) start=atol(s);
			continue;
			}
		else{
			if(fg==0){									//======== variableStep
				if((s=strtok(s," \t\n\r"))==0) continue;
				beg=atoi(s); end=beg+span;
				if((s=strtok(0," \t\n\r"))==0) continue;
				score=atof(s);
			}
			else if(fg==1){										//======== fixedStep
				beg=start; end=beg+span; start=beg+step;
				score=atof(strtok(s," \t\n\r"));
			}
			if(score>0) addSegment(chrom,beg,end,score);
		}
	}
}

void writeProfile(){
	FILE *f=fopen("x","wt");
	for(int i=-tstWindow; i<tstWindow; i++){
		int pp=i+tstWindow;
		fprintf(f,"%i\t%f\n",i,dens[pp]);
	}
	fclose(f);
}

void test1(){
	readPoints("B_G.mRNA_gene_beg_act.bed");
	getMem(dens,tstWindow*2, "test1");
	readWig("F_Br.H3K4me1.wig");
	writeProfile();

//	testFind("chr14",35179587);
//	testFind("chr14",35100000);
//	testFind("chr14",35170000);
//
//
}

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

const char* gtypes[]={
		"pseudogene",
		"lincRNA",
		"protein_coding",
		"antisense",
		"processed_transcript",
		"snRNA",
		"sense_intronic",
		"miRNA",
		"misc_RNA",
		"snoRNA",
		"rRNA",
		"3prime_overlapping_ncrna",
		"polymorphic_pseudogene",
		"sense_overlapping",
		"IG_V_gene",
		"IG_C_gene",
		"IG_J_gene",
		"IG_V_pseudogene",
		"TR_C_gene",
		"TR_J_gene",
		"TR_V_gene",
		"TR_V_pseudogene",
		"IG_C_pseudogene",
		"TR_D_gene",
		"TR_J_pseudogene",
		"IG_J_pseudogene",
		"IG_D_gene",
		"Mt_tRNA",
		"Mt_rRNA",
};
FILE *outf[30];
int nGt=30;

int findGt(char *s){
	for(int i=0; i<nGt; i++){
		if(strcmp(gtypes[i],s)==0) return i;
	}
	return -1;
}

void DecodeGENECODE(){
	char b[4096],attr[4096];
	int i=0;
	verbose=1;
	zeroMem(outf,30);
	FILE *f=fopen("gencode.v19.annotation.gtf","rt");
	for(;fgets(b,sizeof(b),f); ){
		if(b[0]=='#') continue;
		if(i%100000==0) verb("%i\n",i); i++;
		if(i>10000) break;
		char *chr=strtok(b,"\t\n\r");
		char *src=strtok(0,"\t\n\r");
		char *type=strtok(0,"\t\n\r");
		if(strcmp(type,"gene")!=0) continue;
		char *beg=strtok(0,"\t\n\r");
		char *end=strtok(0,"\t\n\r");
		char *score=strtok(0,"\t\n\r");
		char *strand=strtok(0,"\t\n\r");
		char *frame=strtok(0,"\t\n\r");
		char *s=strtok(0,"\t\n\r"); strcpy(attr,s);
		s=strstr(s,"gene_type");
		char *gt=strchr(s,'\"'); if(gt==0) continue;
		gt++;
		s=strchr(gt,'\"'); if(s==0) continue;
		*s=0;
		int k=findGt(gt); if(k<0) continue;
		if(outf[k]==0) {
			char oo[512]; sprintf(oo,"%s.gtf1",gt);
			outf[k]=fopen(oo,"wt");
		}
		fprintf(outf[k],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				chr,src,type,beg,end,score,strand,frame,attr);
	}
	verb("TOTAL=%i\n",i);
	fclose(f);
	for(int i=0; i<30; i++){
		if(outf[i]) fclose(outf[i]);
	}
}

float *tr1,*tr2;

float sumCC(int shift){
	int ll=bTrack1.lProf;
	int i0=0, j0=shift;
	if(j0 <0) {i0-=shift; j0=0;}
	float cc=0;
	for(int i=i0, j=j0; i<ll && j <ll; i++,j++)
		cc+=tr1[i]*tr2[j];
	printf("%i\t%f\n",shift,cc);
	return cc;
}

void checkCrossCorr(){
	int l=bTrack1.lProf;
	int w=200;
	float *cc;
	getMem(cc,w*2,"chkCorr #1"); zeroMem(cc,w*2);
	getMem(tr1,l,"chkCorr #2");
	getMem(tr2,l,"chkCorr #3");
	double av1=0,av2=0,sd1=0,sd2=0;
	for(int i=0; i<l; i++){
		float x;
		x=bTrack1.getValue(i,0);
		tr1[i]=x; av1+=x; sd1+=x*x;
		x=bTrack2.getValue(i,0);
		tr2[i]=x; av2+=x; sd2+=x*x;
	}
	av1/=l; av2/=l;
	sd1=sqrt((sd1-l*av1*av1)/(l-1));
	sd2=sqrt((sd2-l*av2*av2)/(l-1));
	for(int i=0; i<l; i++){
		tr1[i]=(tr1[i]-av1)/sd1;
		tr2[i]=(tr2[i]-av2)/sd2;
	}
	for(int i=-w; i<w; i++){
		cc[i+w]=sumCC(i);
	}
	FILE *f=fopen("CC","wt");
	for(int i=0; i<w*2; i++){
		fprintf(f,"%i\t%f\n", (i-w)*100, cc[i]);
	}
	fclose(f);
}

void test(){
	clearDeb();
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	FILE *xml=0;


	xml=fopen("statistics.xml","r+");
	int ll=fseek(xml,-7,SEEK_END);

	const char*bb="<run> </run> \n</xml>\n";
	fprintf(xml,"<run> </run> \n</xml>\n");

//deb("$$$$$$$$$$$$$$$$$$ %ld",ll);

//fprintf(xml,"<run> </run> \n</xml>\n");
fclose(xml);



//	const char *pcname="proj";
//	const char *b="xxx";
//	double avBg=1,avFg=2.5,sdBg=3.78,sdFg=1.735;
//	statTest *MannW=new statTest();
//	MannW->z=0.35; MannW->pVal=1.e-75;
//	FILE *xml=fopen("t.xml","wt");

	exit(0);
}



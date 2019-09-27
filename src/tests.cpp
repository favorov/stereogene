/*
 * tests.cpp
 *
 *  Created on: 22.02.2013
 *      Author: 1
 */
#include "track_util.h"
//#include <sys/types.h>
//#include <sys/stat.h>
#include <unistd.h>


//struct Histogram{
//	double  minVal, maxVal,  // Min & Max values. Min=-1; Max=1;
//			bin,			 // bin size
//			e,  			 // Mean
//			sigma,			 // standard deviation
//			alpha,beta;		 // parameters fof Beta-distribution
//	int 	nBin,			 // Number of bins
//			count;			 // Number of observations
//	double  *dd,			 // Distribution density
//			*db,			 // Appropriate Beta density
//			*Fp,			 // Cummulative  left distribution  ( F(x) )
//			*Fm;			 // Cummulative  right distribution ( 1-F(x) )
//	int iq,im,iqq;
//
//	bool ready;				 // block add() and norm() after norm()
//	Histogram(int nBin);
//	void add(double x);
//	void norm();
//	void normBeta();		 // Calculate cummulative Beta distribution
//	void normF();			 // Calculate cummulative real distribution
//	double pValp(double x);
//	double pValm(double x);
//	double interpol(double x, double* fun);
//	void print(FILE *f);
//	void calcCDF(double *d);
//};
//
//Histogram::Histogram(int n){
//	minVal=-1; maxVal=1; e=0; sigma=0; alpha=beta=1; iq=iqq=im=0;
//	nBin=n; count=0;
//	dd=0; db=0; Fp=0; Fm=0;
//	errStatus="init Histogram";
//	getMem(dd, nBin, "Histogram #1");
//	getMem(db, nBin, "Histogram #1");
//	getMem(Fp, (nBin+1), "Histogram #1");
//	getMem(Fm, (nBin+1), "Histogram #1");
//	zeroMem(dd,nBin);
//	bin=(maxVal-minVal)/nBin;
//	ready=false;
//	errStatus=0;
//}
//
////=================== Add value to the histogram ===========
//void Histogram::add(double x){
//	if(ready) return;
//	int i=int((x-minVal)/bin);
//	if(i<0) i=0; if(i>=nBin) i=nBin-1;
//	dd[i]++; e+=x; sigma+=x*x;
//	count++;
//}
////============= Calculate cummulative Beta distribution for the background distribution
//void Histogram::normBeta(){
//	if(ready) return; ready=true;
//	norm();
//	double eb=0;
//	for(int i=0; i<nBin; i++){	// calculate the integral of the beta distrib.
//		double x=minVal+bin*i;
//		eb+=(db[i]=xBetaD(alpha,beta,x));
//	}
//	for(int i=0; i<nBin; i++)	{	// normalyze beta distribution
//		db[i]/=eb*bin;
//	}
//	calcCDF(db);
//}
////======================= calculate cumulative distributions
//void Histogram::calcCDF(double *d){
//	Fp[0]=Fm[nBin]=0; Fm[0]=1;
//	for(int i=1, j=nBin-1; i<nBin; i++,j--){
//		Fp[i]=Fp[i-1]+d[i]*bin;
//		Fm[j]=Fm[j+1]+d[j]*bin;
//	}
//
//}
////======================== normalize the foreground distributions
//void Histogram::normF(){
//	if(ready) return; ready=true;
//	norm();
//	calcCDF(dd);
//}
////===================== normalize the histogram and calculate the statistical parameters
//void Histogram::norm(){
//	e/=count; sigma=(sigma-count*e*e)/(count-1);
//	double eBeta=2e-1,  d=sigma/4, gg=(1-eBeta)/eBeta, gg1=gg+1;
//	alpha=(gg/(d*gg1*gg1) - 1)/(gg1);
//	beta=gg*alpha;
//
//	sigma=sqrt(sigma);
//	for(int i=0; i<nBin; i++){
//		dd[i]=dd[i]/count/bin;
//	}
//}
//
////============================================ Interpolation
//double Histogram::interpol(double x, double *fun){
//	int i=int((x-minVal)/bin);
//	if(i<0) {i=0;}
//	if(i>=nBin) {i=nBin-1;}
//	double dx0=x-(minVal+bin*i), dx1=bin-dx0;
//	double f0=log(fun[i]), f1=log(fun[i+1]);
//	return exp((f0*dx1+f1*dx0)/bin);
//}
////=======================================================
//
//double Histogram::pValp(double x){return interpol(x,Fp);}
//double Histogram::pValm(double x){return interpol(x,Fm);}
//
//double xBetaD(double a, double b, double x){
//	return pow(1-x,a)*pow(1+x,b);
//}
//
//void Histogram::print(FILE *f){
//	for(int i=0; i<nBin; i++){
//		double x=minVal+bin*i;
//		fprintf(f,"%f\t%f\t%.2e\t%.2e\t%.2e\n",x,dd[i], db[i],Fp[i],Fm[i]);
//	}
//}
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


double *gauss(int l){
	double *x; getMem(x,l, "gauss");
	double a=0.99;
	x[0]=0;
	for(int i=1; i<l; i++)
		x[i]=a*x[i-1]+rGauss(0,0.1);
	return x;
}

//=========================================================== print chromosome data (testing)
void printChrom(){
   for(int i=0; i<n_chrom; i++){
      Chromosome chrom=chrom_list[i];
      printf("%s\t%li\t%i\n",chrom.chrom,chrom.length,chrom.base);
   }
}

//===========================================================
//void bTrack::printBytes(int from, int to){
//	printBytes(stdout, from,to);
//}
//void bTrack::printBytes(FILE *f, int from, int to){
//	fprintf(f,">>>>>>>>> %i..%i\n",from,to);
//	for(int i=from; i<to; i++) {if(i%25==0) fprintf(f,"\n");fprintf(f," %3i",bytes[i]);}
//	fprintf(f,"\n");
//	if(hasCompl){
//		fprintf(f,"<<<<<<<<<");
//		for(int i=from; i<to; i++)  {if(i%25==0) fprintf(f,"\n");fprintf(f," %3i",cbytes[i]);}
//		fprintf(f,"\n");
//	}
//}
//void bTrack::writeBytes(FILE*f){ // Print bytes (for testing)
//	int n=0, t=15;
//	for(int i=0; i<lProf; i++) {
//		if(bytes[i]>15) {fprintf(f," %3i\n",bytes[i]); n++;}
//	}
//	printf("t=%i n=%i\n",t,n);
//}

//void bTrack::printWindow(int n){
//	printWindow(stdout,n);
//}
//void bTrack::printWindow(FILE *f,int n){
//	fprintf(f,"==========\n");
//	for(int i=0; i<n; i++) fprintf(f," %4.1f",profWindow[i]);
//	fprintf(f,"\n");
//}
//===============================================================================
//===============================================================================
//===============================================================================
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

//===================================================================
//===================================================================
double* GaussProcess(double *w,int  n, int kk){
	getMem(w,n,"tst");
	double a=0.995, b=1;
	double x=0, xmin=1.e+200, xmax=-xmin;
	for(int i=0; i<n; i++){
		x=a*x+b*rGauss();
		double xx=exp(-0.1*x);
		w[i]=xx;
	}

	for(int i=0; i<n; i++){
		if(kk && ((i/kk)%2 == 0)) {
			w[i]+=10;
		}
		double xx=w[i];
		if(xmin > xx) xmin=xx;
		if(xmax < xx) xmax=xx;
	}


	for(int i=0; i<n; i++){
		w[i]=1000*(w[i]-xmin)/(xmax-xmin);
	}
	return w;
}

//==================================================================
int genomeSize=100000;
int lStep_test=1000;
double * random_wig(double *w,int num){
	int genomeSize=100000;

	w=GaussProcess(w,genomeSize, num*lStep_test);

	char b[256]; sprintf(b,"rnd_%02i.wig",num);

	FILE *f=fopen(b,"wt");
	fprintf(f,"track type=wiggle_0 description=\"-sin\"\n");
	fprintf(f,"fixedStep chrom=chr1 start=0 step=100 span=100\n");
		for(int i=0; i<genomeSize; i++){
			int k=(int)w[i];
			fprintf(f,"%i\n",k);
		}
	fclose(f);
	return w;
}

//===========================================================================
int l_chrom=2000000;
int n_nucl = 100000;
int *nucl_pos;
char *nucl_prop;
int bin=30;
int delta=20;
float *pr1, *pr2;
double pre,prd,pr1e,pr2e,pr1d,pr2d;

void makeLabel(float *pr, double lam, int p){
	zeroMem(pr,l_chrom);
	for(int i=0; i<n_nucl; i++){
if(i%1000 ==0) verb("i=%i\r",i);
		int k=nucl_pos[i];
		if(k<0) continue;
		if((nucl_prop[i] & p)==p){
			double x=-log(drand()/lam)*lam;
			double sigma=delta/3.; sigma*=sigma;

			int j0=k-delta, j1=k+delta;
			if(j0<0) j0=0; if(j1>=l_chrom) j1=l_chrom;
			for(int j=j0; j<j1; j++){
				float d=j-k;
				double y=x*exp(-d*d/sigma);
				pr[j]+=y;
			}
		}
	}
	pre=0; prd=0;
	double nz=0;
	for(int i=0; i<l_chrom; i++){
		double x=pr[i];
		if(x==0) nz++;
		pre+=x; prd+=x*x;
	}
	pre/=l_chrom; prd=prd/l_chrom-pre*pre; prd=sqrt(prd);
	verb("e=%f sd=%f   nz=%f\n",pre,prd, nz/l_chrom);
}

void wr(const char *fname, float *pr){
	FILE *f=fopen(fname,"w");
	fprintf(f,"track type=bedGraph name=\"%s\" description=\"test\"\n", fname);
	for(int i=0; i<l_chrom; i++){
if(i%1000000==0) verb("wr %i\r",i);
		ScoredRange pos, pos0;
		double x=int(pr[i]*10)/10.;
		pos.score=x; pos.end=(pos.beg=i*bin)+bin;
		if(pos.score == pos0.score){
			pos0.end=pos.end; continue;
		}
		if(pos0.chrom !=0)
			pos0.printBGraph(f);
		pos0.beg=pos.beg; pos0.end=pos.end; pos0.chrom=(char*)"chr1"; pos0.score=pos.score;
		pos0.printBGraph(f);
	}
verb("\n");
	fclose(f);
}

void wr_chrom(){
	FILE *f=fopen("test_chr","w");
	fprintf(f,"chr1\t%i\n",l_chrom*bin);
	fclose(f);
}

double p0[]={15,40,40, 5};
double pTrans[4][4]={
		{ 3,	47,	47,	3},
		{ 4,	90,	3,	3},
		{ 4,	3,	90,	3},
		{ 2,	2,	2,	94},
};

void normTrans(){
	double xx=0;
	for(int a=0; a<4; a++){
		xx=0;
		for(int i=0; i<4; i++) xx+=pTrans[a][i];
		for(int i=0; i<4; i++) pTrans[a][i]/=xx;
	}
	xx=0;
	for(int i=0; i<4; i++) xx+=p0[i];
	for(int i=0; i<4; i++) p0[i]/=xx;
	printf("---------------------\n");
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++)
			printf("%f\t",pTrans[i][j]);
		printf("\n");
	}
	printf("---------------------\n");
	for(int i=0; i<4; i++) printf("%f\t",p0[i]);
		printf("\n");

}

int next(int a){
	double r=drand(), x=0;

//	for(int i=0; i<4; i++){
//		x+=p0[i];
//		if(r<x) return i;
//	}
//
	x=0;
	for(int i=0; i<4; i++){
		x+=pTrans[a][i];
		if(r < x) return i;
	}
	return 0;
}

int intcmp(const void* a, const void* b){
	return *((int*)a)-*((int*)b);
}

void makeNucleosomes(){
	for(int i=0; i<n_nucl; i++){
		nucl_pos[i]=randInt(l_chrom);
	}
	qsort(nucl_pos,n_nucl,sizeof(int),intcmp);
	for(int i=1; i<n_nucl; i++){
		if(nucl_pos[i]==nucl_pos[i-1]) nucl_pos[i-1]=-1;
	}
	//================================ Markov process
	double x=0, r=drand();
	for(int i=0; i<4; i++){
		x+=r; if(x < p0[i]) nucl_prop[0]=i;
		break;
	}

	int n0=0, n1=0, n2=0, n3=0;

	for(int i=1; i<n_nucl; i++){
		nucl_prop[i]=next(nucl_prop[i-1]);
		int p=nucl_prop[i];
		if(p==0) n0++;
		if(p==1) n1++;
		if(p==2) n2++;
		if(p==3) n3++;
	}

deb(">>>>>>>>> n1=%i  n2=%i  n3=%i",n1,n2,n3);
double nn=n0+n1+n2+n3;
printf("%f\t%f\t%f\t%f\n",n0/nn,n1/nn, n2/nn, n3/nn);

}

void test(){
	char b[256];
	clearDeb();
	debugFg=DEBUG_LOG|DEBUG_PRINT;
verb("================ GENERATE TEST DATA =============\n");
	getMem(pr1,l_chrom,0);
	getMem(pr2,l_chrom,0);
	getMem(nucl_pos,n_nucl,0);
	getMem(nucl_prop,n_nucl,0);
	normTrans();
	wr_chrom();

	makeNucleosomes();
	makeLabel(pr1,8,1); pr1e=pre; pr1d=prd;
	outFile=strcpy(b,"label1");
	wr("l1.bgraph",pr1);

//	makeNucleosomes();
	makeLabel(pr2,8,2); pr2e=pre; pr2d=prd;
	outFile=strcpy(b,"label2");
	wr("l2.bgraph",pr2);
	double cc=0;
	for(int i=0; i<l_chrom; i++){
		cc+=(pr1[i]-pr1e)*(pr2[i]-pr2e)/pr1d/pr2d;
	}
	cc/=l_chrom;
	deb("=======>>>> cc=%f",cc);
	exit(0);
}
//===================================================================
//===================================================================
//===================================================================
//===================================================================
const int progType=0;
void printMiniHelp(){;}
void printProgDescr(){;}


int main(int argc, char **argv) {
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	clearDeb();
	verbose=1;
	return 0;

//	test();
}

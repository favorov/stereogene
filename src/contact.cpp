/*
 * contact.cpp
 *
 *  Created on: Oct 26, 2017
 *      Author: mironov
 */
#include "track_util.h"
#include <unistd.h>

const int progType=0;
void printMiniHelp(){;}
void printProgDescr(){;}
void coverage(const char* fname);
char *qswe;
int aa=0;

char * strtab(char *b){
	char *s, *ss;
	if(b){qswe=b;}
	ss=s=qswe;
	for(;*ss; ss++){
		if(*ss == '\t' || *ss=='\n'){
			*ss=0; qswe=ss+1;
			return s;
		}
	}
	return s;
}

struct VOCAB_ENTRY{
	char *name;
	int count;
	VOCAB_ENTRY(const char *nm);
};

struct GENE_ENTRY:VOCAB_ENTRY{
	VOCAB_ENTRY *type;
};

struct VOCAB{
	VOCAB_ENTRY **ent;
	int n;
	int capacity;
	VOCAB();
	VOCAB_ENTRY*get(const char* nm);
	VOCAB_ENTRY*add(const char* nm);
	void sort();
	void db(){
		deb("==============================================================  n=%i",n);
		for(int i=0; i<n; i++) deb("-> %6i %s", ent[i]->count, ent[i]->name);
	}
};

VOCAB *geneNames;
VOCAB *geneTypes;
VOCAB *chrom;

struct DNA{
	char *chr;
	int pos;
	int fg;
};

struct LOCUS{
	VOCAB_ENTRY *chr;
	int b,e;
	char strand;
	void read();
};

void LOCUS::read(){
	chr=chrom->add(strtab(0));
	b=atoi(strtab(0));
	e=atoi(strtab(0));
	strand=*((strtab(0)));
}

VOCAB_ENTRY::VOCAB_ENTRY(const char *nm){
	name=strdup(nm);
	count=0;
}

VOCAB::VOCAB(){
	capacity=40000;
	getMem(ent,capacity,0);
	n=0;
}

int ent_cmp(const void *e1, const void *e2){
	VOCAB_ENTRY **ve1=(VOCAB_ENTRY **)e1, **ve2=(VOCAB_ENTRY **)e2;
	return (*ve2)->count - (*ve1)->count;
}

void VOCAB::sort(){
	qsort(ent, n, sizeof(VOCAB_ENTRY*), ent_cmp);
}

VOCAB_ENTRY*VOCAB::get(const char* nm){
	for(int i=0; i<n; i++){
		if(strcmp(ent[i]->name,nm)==0) return ent[i];
	}
	return 0;
}

VOCAB_ENTRY*VOCAB::add(const char* nm){
	if(*nm==0) nm="_";
	VOCAB_ENTRY *e=get(nm);
	if(e==0){
		e=new VOCAB_ENTRY(nm);
		if(n>=capacity-1){
			capacity+=10000;
			realocMem(ent,capacity,0);
		}
		ent[n++]=e;
	}
	e->count++;
	return e;
}

//number
//rna1_chr rna1_start rna1_end rna1_strand
//rna1_gene_strand rna1_gene_type rna1_gene_name
//dna_chr dna_start dna_end dna_gene_s

struct CONTACT{
	int num;
	LOCUS r1;
	char *rna_gene_name;
	char *rna_gene_type;
	LOCUS dna;
	void read(char *b);
	CONTACT(char *b){read(b);}
};

void CONTACT::read(char *b){
	char *s=strtab(b);
	num=atoi(s);
	r1.read();
	strtab(0); //skip rna1_gene_strand
	rna_gene_type=geneTypes->add(strtab(0))->name;
	rna_gene_name=geneNames->add(strtab(0))->name;
	dna.read();
}

const int MAX_CONT=50000000;
CONTACT *contacts[MAX_CONT];
int n_cont=0;

const char *tt[]={
		"_",
		"NA",
		"MALAT1",
		"RNU2-2P",
		"LINC00534",
		"RP11-6N13.1",
		"RP11-328J2.1",
		"NEAT1",
		"RP11-124N3.3",
		"RP11-727F15.12",
		"RP11-727F15.13",
		"RNU6-118P",
		"RN7SL119P",
		"RP11-727F15.11",
		"CTD-2374C24.1",
		"AC104389.28",
		"AL109763.2",
		"AC096559.1",
		"GAS5",
		"MBNL1",
		"SNHG1",
		"C22orf34",
		"RAD51B",
		"RNU12",
		"RP11-317N12.1",
		"AC002539.1",
		"CCDC26",
		"DLEU1",
		0
};
const char *gt[]={
//		"3prime_overlapping_ncrna",
//		"antisense",
//		"IG_C_pseudogene",
//		"IG_V_gene",
//		"IG_V_pseudogene",
//		"lincRNA",
//		"miRNA",
//		"misc_RNA",
//		"Mt_rRNA",
//		"Mt_tRNA",
//		"polymorphic_pseudogene",
//		"processed_transcript",
//		"protein_coding",
//		"pseudogene",
//		"rRNA",
//		"sense_intronic",
//		"sense_overlapping",
		"snoRNA",
		"snRNA",
//		"TR_J_gene",
		0
};
//========================================= read contacts ==========================
void writeContacts(const char *name, int fg){
	char b[4096];
	if(*name==0) name="_";
	sprintf(b,"data/%s_Gav2_AD08.bed",name);
	FILE * out=fopen(b,"w");
	verb("write %s ...\n",b);

	for(int i=0; i<n_cont; i++){
		CONTACT *c=contacts[i];
		char *nm=fg ? c->rna_gene_name : c->rna_gene_type;
		if(strcmp(nm,name)==0){
			fprintf(out,"%s\t%i\t%i\t1\n",c->dna.chr->name,c->dna.b,c->dna.e);
		}
	}
	fclose(out);
}
void readContacts(const char *fname){
	geneTypes=new VOCAB();
	geneNames=new VOCAB();
	chrom=new VOCAB();
	FILE *f=xopen(fname,"r");
	char b[4096];
	fgets(b,sizeof(b),f);	//skip header
	for(int i=0; fgets(b,sizeof(b),f); i++){
		if(i%100000==0) verb("%i\n",i/1000);
//if(i>1000000) break;
		contacts[n_cont++]=new CONTACT(b);
		if(n_cont >= MAX_CONT) {
			verb("too many contacts\n");
			exit(1);
		}
	}
	geneTypes->sort();
	geneNames->sort();

//	for(int i=0; tt[i]; i++){
//		writeContacts(tt[i],1);
//	}
	for(int i=0; gt[i]; i++){
		writeContacts(gt[i],0);
	}
}


//====================================== MAKE COVERAGE ====================
int dn_cmp(const void* d1, const void* d2){
	DNA *dd1=(DNA*)d1, *dd2=(DNA*)d2;
	int a=strcmp(dd1->chr, dd2->chr);
	if(a) return a;
	a=dd1->pos-dd2->pos; if(a) return a;
	return dd1->fg-dd2->fg;
}

void coverage(const char* fname){
	char b[4096], *s;

//	deb(getcwd(b,sizeof(b)));
	verb("%s\n", fname);
	sprintf(b,"data/%s_Gav2_AD08.bed", fname); FILE *in =xopen(b,"r");
	sprintf(b,"data/%s_Gav2_AD08.bgr", fname); FILE *out=xopen(b,"w");
	sprintf(b,"%s.ll" , fname); FILE *ft =xopen(b,"w");

	DNA *dn;
	int n_dn=0;
	getMem(dn,60000000,0);

	for(int i=0;(s=fgets(b,sizeof(b),in))!=0; i++){
		char *chr;
		int beg, end;
		if(i%100000==0) verb("%i\n",i/1000);
		char bb[256]; strcpy(bb,strtok(b,"\n"));
		chr=strtok(b,"\t\n"); 			//chrom
		beg=atoi(strtok(0,"\t\n"));		// dna start
		end=atoi(strtok(0,"\t\n"));		// dna end
		chr=strdup(chr);
		dn[n_dn].chr=chr; dn[n_dn].pos=beg; dn[n_dn].fg=1; n_dn++;
		dn[n_dn].chr=chr; dn[n_dn].pos=end; dn[n_dn].fg=0; n_dn++;
	}
	verb("sort...\n");
	qsort(dn,n_dn,sizeof(DNA),dn_cmp);

	verb("write...\n");
	fprintf(out,"track type=bedGraph name=\"%s\" description=\"%s\"\n",fname,fname);

	int count=0, pos=0;
	for(int i=0 ;i < n_dn; i++){
				if(count && pos!=dn[i].pos){
					fprintf(out,"%s\t%i\t%i\t%i\n",dn[i].chr, pos,dn[i].pos,count);
					if(count >= 50)
						fprintf(ft,"%s\t%i\t%i\n",dn[i].chr,pos, count);
				}
				if(dn[i].fg) {count++; pos=dn[i].pos;}
				else 		 {count--; pos=dn[i].pos;}
	}
	fclose(in); fclose(out); fclose(ft);
	free(dn);
	verb("done\n");
}
//========================================================
int getNameIdx(char *nm, char **names, int &n){
	char b[4096];
	nm=strcpy(b,nm);
	if(memcmp(nm,"confounder.proj/",5)==0) nm+=16;
	if(memcmp(nm,"hg19_",5)==0) nm+=5;
	if(memcmp(nm,"K562_",5)==0) nm+=5;
	char *s=strstr(nm,"_Gav2_AD08.bgr"); if(s) *s=0;
	s=strstr(nm,"_100.lifted.bed"); if(s) *s=0;
	s=strstr(nm,".bgr"); if(s) *s=0;
	s=strstr(nm,".broadPeak"); if(s) *s=0;

	for(int i=0; i<n; i++){
		if(strcmp(names[i],nm)==0) return i;
	}
	names[n++]=strdup(nm);
	return n-1;
}


void readStat(){
	FILE *f=fopen("statistics","r");
	char *hg19[500];
	char *gav2[500];
	double mtx[500][500];
	int n_hg19=0, n_gav2=0;
	char b[4096];
//Fg_av_Corr	FgCorr_sd	Bg_Corr 	Bg_av_Corr	BgCorr_sd	Mann-Z  	p-value	pcorProfile
	fgets(b,sizeof(b),f);
	for(;fgets(b,sizeof(b),f);){
		char *s;
		int i=0,j=0;
		s=strtok(b,"\t\n"); //skip id
		s=strtok(0,"\t\n"); //skip Date
		s=strtok(0,"\t\n"); //skip version
		s=strtok(0,"\t\n"); //get name1
		//==========================================
		if(strstr(s,"hg19_")) i=getNameIdx(s,hg19, n_hg19);
		else 				  j=getNameIdx(s,gav2, n_gav2);
		//==========================================
		s=strtok(0,"\t\n"); //get name2
		if(strstr(s,"hg19_")) i=getNameIdx(s,hg19, n_hg19);
		else 				  j=getNameIdx(s,gav2, n_gav2);
		//==========================================
		s=strtok(0,"\t\n"); //skip wSize
		s=strtok(0,"\t\n"); //skip kernelType
		s=strtok(0,"\t\n"); //skip nFgr
		s=strtok(0,"\t\n"); //skip nBkg
		s=strtok(0,"\t\n"); //get Fg_Corr
		double corr=atof(s);
		mtx[i][j]=corr;
	}
	fclose(f);
	f=fopen("tb","w");
	fprintf(f,"k562\t");
	for(int j=0; j<n_gav2; j++){
		fprintf(f,"%s\t",gav2[j]);
	}
	fprintf(f,"\n");
	for(int i=0; i<n_hg19; i++){
		fprintf(f,"%s\t",hg19[i]);
		for(int j=0; j<n_gav2; j++){
			fprintf(f,"%0.3f\t",mtx[i][j]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
	exit(0);
}

int main(int argc, char **argv) {
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	clearDeb();
	verbose=1;

	readStat();
//	readContacts("Gav2_AD08.all.contacts.tab");
//	exit(0);

	for(int i=0; gt[i]; i++){
		coverage(gt[i]);
	}

//	for(int i=0;tt[i];i++){
//		coverage(tt[i]);
//	}

	return 0;
}






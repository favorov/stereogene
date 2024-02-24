#include "track_util.h"


//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <cmath>
//#include <ctype.h>


char inFile[4096];
FILE *fin;
FILE *fgene  ;
FILE *fg_beg ;
FILE *fg_end ;
FILE *fexn   ;
FILE *fe_beg ;
FILE *fe_end ;
FILE *fivs   ;
FILE *fi_beg ;
FILE *fi_end ;




const int progType=PG;
const char * progName="parse_genes";






struct nameList{
	char ** names;
	int nNames;
	int maxNames;
	char *getName(char *s);
	nameList();
};


nameList::nameList(){
	maxNames=100;
	names=(char **)malloc(maxNames*sizeof(char*));
	nNames=0;
}


char * nameList::getName(char *s){
	for(int i=0; i<nNames; i++){
		if(strcmp(names[i],s)==0) return names[i];
	}
	if(nNames >= maxNames-2){
		maxNames+=1000;
		names=(char **)realloc(names,maxNames*sizeof(char*));
	}
	return (names[nNames++]=strdup(s));
}


nameList chromosomes=nameList();
nameList genetypes  =nameList();




long readInt_pg(){
	char *s=strtok(0," \t\n\r");
	if(s==0) return -1;
	if(! isdigit(*s)) return -1;
	return atol(s);
}
double readFloat_pg(){
	char *s=strtok(0," \t\n\r");
	if(s==0) return -1;
	return atof(s);
}


//================================ PARSE REFSEQ.BED ==========================
//============================================================================
struct pg_refseqIV{
	char *chr;
	int f,t;
	char strand;
	char *name;
	void print(FILE *f);
};


void pg_refseqIV::print(FILE *ff){
	fprintf(ff,"%s\t%i\t%i\t%s\t1\t%c\n",chr,f,t, name ,strand);
}




const int CAPAS=100000;
struct pg_refseqSet{
	char *setName;
	pg_refseqIV *ivs;
	int n;
	int capacity;
	pg_refseqSet(const char *sn);
	void addIV(char *chrom, int from, int to, const char strand, const char *name);
	void print(FILE*f);
};


pg_refseqSet::pg_refseqSet(const char *sn){
	setName=strdup(sn);
	capacity=CAPAS;
	ivs=(pg_refseqIV*) malloc(capacity*sizeof(pg_refseqIV));
	n=0;
}


void pg_refseqSet::addIV(char *chrom, int from, int to, const char strand, const char *name){


	if(n > 0){
		pg_refseqIV *prev=ivs+n-1;
		if(prev->chr == chrom && prev -> f <= to && prev->t >= from && prev->strand==strand){
			if(prev->f >= from) prev->f=from;
			if(prev->t <= to  ) prev->t=to;
			return;
		}
	}
	if(n >= capacity-2){
		capacity+=CAPAS;
		ivs=(pg_refseqIV*) realloc(ivs,capacity*sizeof(pg_refseqIV));
	}
	pg_refseqIV *iv=ivs+n;
	iv->chr=chrom;
	iv->f=from;
	iv->t=to;
	iv->strand=strand;
	iv->name=strdup(name);
	n++;
}


int rsIvCmp(const void *r1, const void *r2){
	pg_refseqIV *iv1=(pg_refseqIV *) r1;
	pg_refseqIV *iv2=(pg_refseqIV *) r2;


	int x=strcmp(iv1->chr, iv2->chr);
	if(x) return x;
	if(iv1->f == iv2->f) {
		if(iv1->strand != iv2->strand) return iv1->strand - iv2->strand;
		return iv1->t-iv2->t;
	}
	return iv1->f-iv2->f;
}


void pg_refseqSet::print(FILE *f){
	qsort(ivs,n,sizeof(pg_refseqIV),rsIvCmp);
	for(int i=0; i<n-1; i++){
		pg_refseqIV *cur=ivs+i;
		do{
			pg_refseqIV *next  =ivs+i+1;
			if(cur->chr == next->chr && cur -> f <= next->t && cur->t >= next->f && cur->strand==next->strand){
				if(cur->f >= next->f ) cur->f=next->f;
				if(cur->t <= next->t ) cur->t=next->t;
				i++;
				if(i>=n) break;
			}
			else{
				break;
			}
		}while(true);
		cur->print(f);
	}
}




pg_refseqSet rs_Genes=pg_refseqSet("genes");
pg_refseqSet rs_G_beg=pg_refseqSet("g_beg");
pg_refseqSet rs_G_end=pg_refseqSet("g_end");
pg_refseqSet rs_Exn  =pg_refseqSet("exn");
pg_refseqSet rs_Exn_b=pg_refseqSet("e_beg");
pg_refseqSet rs_Exn_e=pg_refseqSet("e_end");
pg_refseqSet rs_Ivs  =pg_refseqSet("ivs");
pg_refseqSet rs_Ivs_b=pg_refseqSet("i_beg");
pg_refseqSet rs_Ivs_e=pg_refseqSet("i_end");




void parseRefSeq(){
	int SG_BUFEXT=0x100000;
	char *bb=(char*)malloc(SG_BUFEXT), *name, *sx, strand;
	char *inputString;
	for(int iLine=0; (inputString=fgets(bb,SG_BUFEXT,fin))!=0; iLine++){
		if(iLine % 10000 == 0) verb ("%i\n", iLine);
		if(strncmp(inputString,"browser"  ,7)==0) continue;
		if(strncmp(inputString,"track"    ,5)==0) continue;
		if(*inputString=='#' || *inputString==0) continue;


		char *chrom=chromosomes.getName(strtok(inputString,"\t\n"));
		int beg=readInt_pg();									//======== find beg field
		int end=readInt_pg();
		name=strtok(0,"\t\n");
		if(name==0) {
			rs_Genes.addIV(chrom,beg,end,'.',".");
			rs_G_beg.addIV(chrom,beg,beg+1,'.',".");
			rs_G_end.addIV(chrom,end,end+1,'.',".");
			continue;
		}


		sx=strtok(0,"\t\n"); if(sx==0) continue;			//======== take score field
		sx=strtok(0,"\t\n"); if(sx==0) continue;			//======== take strand field
		strand = (*sx==0) ? '.' : *sx;
		rs_Genes.addIV(chrom,beg,end,strand,name);
		if(strand=='-'){
			rs_G_end.addIV(chrom,beg,beg+1,strand,name);
			rs_G_beg.addIV(chrom,end,end+1,strand,name);
		}
		else{
			rs_G_beg.addIV(chrom,beg,beg+1,strand,name);
			rs_G_end.addIV(chrom,end,end+1,strand,name);
		}
//====================================== genes ================================
		sx=strtok(0,"\t"); if(sx==0) continue;			//======== skip tick1
		sx=strtok(0,"\t"); if(sx==0) continue;			//======== skip tick2
		sx=strtok(0,"\t"); if(sx==0) continue;			//======== skip rgb
		sx=strtok(0,"\t"); if(sx==0) continue;			//======== number of exons


		//======================= read gene structure
		int nn=atoi(sx);
		if(nn<2) {
			rs_Exn.addIV(chrom, beg, end, strand, name);
			rs_Exn_b.addIV(chrom, end, end+1, strand, name);
			rs_Exn_e.addIV(chrom, beg, beg+1, strand, name);
			continue;
		}
		int lExn[nn];
		int posExn[nn];


		for(int i=0; i<nn; i++) {
			sx=strtok(0,",\t");
			if(sx==0){beg=-1; break;}
			lExn[i]=atoi(sx);
		}


		for(int i=0; i<nn; i++) {
			sx=strtok(0,",\t");
			if(sx==0){beg=-1; break;}
			posExn[i]=atoi(sx);
		}


		int genePos=beg;
//===== now 'beg' and 'end' are positions relative to gene beg/
//===== exons ================================================
		for(int i=0; i<nn; i++){
			end=(beg=posExn[i])+lExn[i];
			beg=genePos+beg; end=genePos+end;
			rs_Exn.addIV(chrom, beg, end, strand, name);
			if(strand=='-'){
				rs_Exn_e.addIV(chrom, beg, beg+1, strand, name);
				rs_Exn_b.addIV(chrom, end, end+1, strand, name);
			}
			else{
				rs_Exn_b.addIV(chrom, beg, beg+1, strand, name);
				rs_Exn_e.addIV(chrom, end, end+1, strand, name);
			}
		}
//===== introns ==============================================
		for(int i=0; i<nn-1; i++){
			beg=posExn[i]+lExn[i];
			end=posExn[i+1];
			beg=genePos+beg; end=genePos+end;
			rs_Ivs.addIV(chrom, beg, end, strand, name);
			if(strand=='-'){
				rs_Ivs_e.addIV(chrom, beg, beg+1, strand, name);
				rs_Ivs_b.addIV(chrom, end, end+1, strand, name);
			}
			else{
				rs_Ivs_b.addIV(chrom, beg, beg+1, strand, name);
				rs_Ivs_e.addIV(chrom, end, end+1, strand, name);
			}
		}
	}






	verb("print genes\n"); rs_Genes.print(fgene);
	verb("print g_beg\n"); rs_G_beg.print(fg_beg);
	verb("print g_end\n"); rs_G_end.print(fg_end);
	verb("print exons\n"); rs_Exn.print(fexn);
	verb("print e_beg\n"); rs_Exn_b.print(fe_beg);
	verb("print e_end\n"); rs_Exn_e.print(fe_end);
	verb("print introns\n"); rs_Ivs.print(fivs);
	verb("print i_beg\n"); rs_Ivs_b.print(fi_beg);
	verb("print i_end\n"); rs_Ivs_e.print(fi_end);
}
//============================================================================
const int PG_TRASCR=1;
const int PG_EXON=3;


struct GenecodeRecord{
	char *chrom;
	int type;
	int from,to;
	char strand;
//	char *transcrID;
	char geneName[1024];
	char *genetype;
	int level;
	int exonNo;
	int parse(char *b);


	void print(FILE *f);
};


void GenecodeRecord::print(FILE *f){
	fprintf(f,"%s\t%i\t%i\t%i\t%c\t%s\t%s\t%i\t%i\n",
			chrom,type,from,to,strand,geneName,genetype,level,exonNo);
}




int GenecodeRecord::parse(char *b){
	char *sx;
	exonNo=-1; level=-1;
	sx=strtok(b,"\t\n");
	chrom=chromosomes.getName(sx);
	strtok(0,"\t\n");		// skip Havana
	sx=strtok(0,"\t\n");	// read type
	if(strcmp(sx,"transcript")==0) type=PG_TRASCR;
	else if(strcmp(sx,"exon")==0) type=PG_EXON;
	else return 0;
	sx=strtok(0,"\t");	from=atoi(sx); // read from
	sx=strtok(0,"\t");	to  =atoi(sx); // read to
	strtok(0,"\t");		// skip '.'
	sx=strtok(0,"\t");	strand=*sx;	// strand
	strtok(0,"\t");		// skip '.'


	sx=strtok(0,"\t");		// read description
	sx=strtok(sx,";");


	for(;sx!=0; sx=strtok(0,";")){
		while(*sx==' ') sx++;
		char *name;
		char *value;
		char *s=strchr(sx,' ');
		if(s)*s=0;
		else continue;
		name=trim(sx);
		value=trim(s+1);
		//
		if(strcmp(name,"gene_type"    )==0) genetype =genetypes. getName(value);
		if(strcmp(name,"gene_name"    )==0) strcpy(geneName,value);
		if(strcmp(name,"exon_number"  )==0) exonNo=atoi(value);
		if(strcmp(name,"level"        )==0) level=atoi(value);


	}
	return 1;
}


struct gp_interval{
	int f,t;
	char strand;


};
struct gp_intervals{
	int begs[0x100000];
	int ends[0x100000];
	gp_interval interv[0x10000];
	int n_beg;
	int n_end;
	int n_ivs;
	void addIv(int f, int t, char strand);
	gp_intervals(){clean();}
	void gp_addPos(int v, int pos[], int &n);
	void clean(){n_beg=0; n_end=0; n_ivs=0;}
	void print(FILE *fiv, FILE *fbeg,  FILE *fend, char *chrom, char *name, char strand);
	void printPos(FILE *f,  int pos[], int n, char *chrom, char *name, char strand);
};


void gp_intervals::printPos(FILE *f, int pos[], int n, char *chrom, char *name, char strand){
	for(int i=0; i<n; i++){
		fprintf(f,"%s\t%i\t%i\t%s\t1\t%c\n",chrom,pos[i],pos[i]+1, name ,strand);
	}
}


int pg_iv_cmp(const void* ivx1, const void *ivx2){
	gp_interval* iv1=(gp_interval*)ivx1;
	gp_interval* iv2=(gp_interval*)ivx2;
	if(iv1->f == iv2->f) return iv1->t - iv2->t;
	return iv1->f - iv2->f;
}
int pg_pos_cmp(const void* px1, const void *px2){
	int *p1=(int*)px1;
	int *p2=(int*)px2;
	return *p1-*p2;
}


void gp_intervals::print(FILE *fiv, FILE *fbeg,  FILE *fend, char *chrom, char *name, char strand){
	qsort(interv, n_ivs,sizeof(gp_interval),pg_iv_cmp);
	qsort(begs, n_beg,sizeof(int),pg_pos_cmp);
	qsort(ends, n_end,sizeof(int),pg_pos_cmp);
	for(int i=n_ivs-1; i > 0; i--){
		if(interv[i].f <= interv[i-1].t  && interv[i].t >= interv[i-1].f){
			if(interv[i].t > interv[i-1].t)
				interv[i-1].t=interv[i].t;
			for(int j=i; j < n_ivs-1; j++){
				interv[j].f=interv[j+1].f;
				interv[j].t=interv[j+1].t;
			}
			n_ivs--;
		}
	}
	for(int i=0; i<n_ivs; i++){
		fprintf(fiv,"%s\t%i\t%i\t%s\t1\t%c\n",chrom,interv[i].f,interv[i].t, name,strand);
	}
	printPos(fend,ends, n_end, chrom, name, strand);
	printPos(fbeg,begs, n_beg, chrom, name, strand);
}






void gp_intervals::gp_addPos(int v, int pos[], int &n){
	for(int i=0; i<n; i++){
		if(pos[i]==v) return;
	}
	pos[n]=v; n++;
}
void gp_intervals::addIv(int f, int t, char strand){
	if(f > t){int x=f; f=t; t=x;}
	if(strand == '-') {gp_addPos(t,begs, n_beg); gp_addPos(f,ends, n_end);}
	else 			  {gp_addPos(f,begs, n_beg); gp_addPos(t,ends, n_end);}
	bool fg=false;
	for(int i=0; i<n_ivs; i++){
		if(interv[i].f <= t && interv[i].t >= f){
			if(interv[i].f >= f) interv[i].f = f;
			if(interv[i].t <= t) interv[i].t = t;
			fg=true; break;
		}
	}
	if(!fg){
		interv[n_ivs].f = f;
		interv[n_ivs].t = t;
		interv[n_ivs].strand=strand;
		n_ivs++;
	}
}


gp_intervals givTranscr=gp_intervals();
gp_intervals givExn    =gp_intervals();
gp_intervals givIvs    =gp_intervals();
void parseGTF(){
	int SG_BUFEXT=0x100000;
	char *bb=(char*)malloc(SG_BUFEXT);
	char *inputString;
	GenecodeRecord record;


	char curGName[256]="";
	char curChr[256]="";
	char curStrand='.';
	int ne=0;
	int exn_end=0;


	for(int iLine=0; (inputString=fgets(bb,SG_BUFEXT,fin))!=0; iLine++){
		if(iLine % 10000 ==0) printf("%i\n", iLine);
		inputString=trim(inputString);
		if(*inputString=='#') continue;


		if(record.parse(inputString)){
//			record.print(fgene);
			if(record.level < pgLevel) {continue;}
			if(strcmp(curGName, record.geneName)){
				givTranscr.print(fgene, fg_beg, fg_end, curChr, curGName, curStrand);
				givExn.print(fexn, fe_beg, fe_end, curChr, curGName, curStrand);
				givIvs.print(fivs, fi_beg, fi_end, curChr, curGName, curStrand);




				givTranscr.clean(); givExn.clean(); givIvs.clean();
				strcpy(curGName,record.geneName);
				strcpy(curChr, record.chrom);
				curStrand=record.strand;
			}
			if(record.type == PG_TRASCR){
				//====================================== make introns
				givTranscr.addIv(record.from, record.to, record.strand);
				ne=0;
			}
			if(record.type == PG_EXON){
				givExn.addIv(record.from, record.to, record.strand);
				if(ne){
					givIvs.addIv(exn_end,record.from, record.strand);
				}
				ne++; exn_end=record.to;
			}
		}
	}
}


//============================================================================
//============================================ Print Help page =========================================
void printMiniHelp(){
	printf("\n");
	printf("The parse_genes program creates bed files for genes/exons/introns starts/bodies/ends\n");
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./parse_genes [-parameters] [RefSeq or GENECODE file]\n");
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}


void printProgDescr(){
	printf("\n");
	printf("The parse_genes program creates bed files for genes/exons/introns starts/bodies/ends\n");
	printf("Usage:\n");
	printf("$ ./parse_genes [-parameters] [RefSeq or GENECODE file]\n");
	printf("\n");
}




int main(int argc, char **argv) {
	char b[5100];
	verbose=true;
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	parseArgs(argc, argv);
	if(nfiles < 1) {printf("no input file defined"); return 0;}
	strcpy(inFile, trim(files[0].fname));


	fin=xopen(argv[1], "r");
	int type=0;	//======================== gencode
	char *s = strrchr(inFile,'.');
	if(s) {
		*s=0; if(keyCmp(s+1,"bed")==0) type=1;		//================ refseq
	}
	sprintf(b,"%s_gene.bed" ,inFile);  fgene =fopen(b,"w");
	sprintf(b,"%s_g_beg.bed",inFile);  fg_beg=fopen(b,"w");
	sprintf(b,"%s_g_end.bed",inFile);  fg_end=fopen(b,"w");
	sprintf(b,"%s_exn.bed"  ,inFile);  fexn  =fopen(b,"w");
	sprintf(b,"%s_e_beg.bed",inFile);  fe_beg=fopen(b,"w");
	sprintf(b,"%s_e_end.bed",inFile);  fe_end=fopen(b,"w");
	sprintf(b,"%s_ivs.bed"  ,inFile);  fivs  =fopen(b,"w");
	sprintf(b,"%s_i_beg.bed",inFile);  fi_beg=fopen(b,"w");
	sprintf(b,"%s_i_end.bed",inFile);  fi_end=fopen(b,"w");
	if(type) parseRefSeq();
	else parseGTF();


	verb("close files\n");
	fclose(fin);
	fclose(fgene);
	fclose(fg_beg);
	fclose(fg_end);
	fclose(fexn);
	fclose(fe_beg);
	fclose(fe_end);
	fclose(fivs);
	fclose(fi_beg);
	fclose(fi_end);
}



/*
 * file_list.cpp
 *
 *  Created on: Apr 3, 2013
 *      Author: mironov
 */

#include "track_util.h"
#define TRACK_BUFF_SIZE 40960


class TrackFile{
public:
	FILE* file;
	char *name,*fname;
	unsigned char *buff;
	int buffsize;
	long curpos;
	long pLength;
	int step;
	int hasCompl;
	int pPcaSgm;
	float* prof;
	TrackFile();
	TrackFile(char *fname);
	float getProfile(int pos);
	void printH();
	int getFilePos(unsigned long pos);
};
int nFiles=0;
TrackFile **tracks;

unsigned long *positions;

long profLength=0;
FILE *outPCA;


TrackFile::TrackFile(){
	file=0; buff=0; buffsize=TRACK_BUFF_SIZE; curpos=0; pLength=0; fname=name=0;step=0;hasCompl=0;prof=0;pPcaSgm=100;
}
void TrackFile::printH(){
	fprintf(outPCA,"\t%s",alTable.convert(name));
}
int TrackFile::getFilePos(unsigned long pos){ return (int)((double)pos*stepSize/step);}

TrackFile::TrackFile(char *fname){
	buffsize=TRACK_BUFF_SIZE;
	errStatus="init Track file";
	getMem0(buff,buffsize, "init Track file #1");

	char bx[4096], bp[4096], b[1024],*s0,*s1;
	strcpy(bx,fname); s0=strrchr(bx,'/'); if(s0==0) s0=bx;
	s1=strrchr(s0,'.'); if(s1) *s1=0;

	strcat(strcpy(bx,profPath),fname);
	strcpy(bp,bx); s0=strrchr(bp,'/'); if(s0==0) s0=bp;
	s1=strrchr(s0,'.'); if(s1) *s1=0;
	strcat(bp,".prm");
	file=xopen(bp,"rt");
	name=0;
	for(;fgets(b,sizeof(b),file)!=0;){
		char *s1=strtok(b,"=");
		char *s2=strtok(0,"=\r\n");
		if(strcmp(s1,"trackName")==0 && s2!=0 && strlen(s2)!=0) name=strdup(skipSpace(s2));
		else if(strcmp(s1,"step"   )==0) step=atoi(s2);
		else if(strcmp(s1,"strand" )==0) hasCompl=atoi(s2);
	}
	fclose(file);
	if(name==0 || strlen(name)==0)
	{
		char b[1024];
		getFnameWithoutExt(b,fname);
		name=strdup(b);
	}
	else 	name=strtok(name," \t");

	name=alTable.convert(name);
	pPcaSgm=pcaSegment/step;
	file=xopen(bx,"rt"); curpos=0;
	fseek(file, 0L, SEEK_END); pLength = (int) ftell(file);
	fseek(file, 0L, SEEK_SET);
	getMem0(prof,nPca, "init Track file #2");
	for(int i=0; i<nPca; i++){
		unsigned long k=(int)( (double)positions[i]*stepSize/step);
		prof[i]=getProfile(k);
	}
	errStatus=0;
}


float TrackFile::getProfile(int pos){
	if(pos < curpos || pos+pPcaSgm > curpos+buffsize){//read buffer
		curpos=(pos/4096)*4096;
		fseek(file,curpos,SEEK_SET);
		fread(buff,1,buffsize,file);
	}
	int pp=pos-curpos;
	float d=0;
	for(int i=0; i<pPcaSgm;i++,pp++){
		if(buff[i]==0) buff[i]=1;
		d+=buff[i]-1;
	}
	d/=pPcaSgm;
	return d;
}

int intCmp(const void*a, const void*b){
	return *(int*)a-*(int*)b;
}

void pcaMain(const char *fname){

	char * names[500];
	char b[1024], *s;

	writeLog("PCA prepare:  %s\n",fname);

	FILE *in=xopen(fname,"rt");
	strcpy(b,fname); s=strrchr(b,'/'); if(s==0) s=b;
	s=strrchr(s,'.'); if(s!=0) *s=0;
	strcat(b,".pca");
	char *outFname=strdup(b);
	outPCA=xopen(b,"wt");

	for(;fgets(b,sizeof(b),in);){
		s=strtok(skipSpace(b),"\n ");
		if(s!=0 && strlen(s)>0 && *s!='#') names[nFiles++]=strdup(s);
	}
	fclose(in);
	//===================== Make positions
	getMem0(positions,nPca, "PCA #1");

	for (int i = 0; i < nPca; i++) {
		positions[i]=randInt(profileLength) ;
	}

	qsort(positions,nPca,sizeof(int),intCmp);

	getMem0(tracks,nFiles, "PCA #2");

	for(int i=0; i<nFiles; i++){
		tracks[i]=new TrackFile(names[i]);
		verb("read profile <%s>\n",tracks[i]->name);
	}

	fprintf(outPCA,"cr\tbeg\tend");
	for(int k=0; k< nFiles; k++){
		tracks[k]->printH();
	}
	fprintf(outPCA,"\n");
	ScoredRange sr;
	int nn=0;
	for(int i=0; i<nPca; i++){
		int k=tracks[0]->getFilePos(positions[i]);
		filePos2Pos(k,&sr,stepSize);

		double z=0;
		for(int k=0; k< nFiles; k++) z=tracks[k]->prof[i];
		if(z!=0){
			nn++;
			fprintf(outPCA,"%s\t%li\t%li",sr.chrom,sr.beg,sr.beg+pcaSegment);
			for(int k=0; k< nFiles; k++){
				fprintf(outPCA,"\t%.1f",tracks[k]->prof[i]);
			}
			fprintf(outPCA,"\n");
		}
	}

	if(RScriptFg){
		strcpy(b,fname); s=strrchr(b,'/'); if(s==0) s=b;
		s=strrchr(s,'.'); if(s!=0) *s=0;
		strcat(b,".r");

		FILE *f=xopen(b,"wt");

		fprintf(f,"aa=read.table(\"%s\",header=TRUE)\n",outFname);
		fprintf(f,"print(cor(aa[,4:%i]))\n",nFiles+3);
		fprintf(f,"pca=prcomp(aa[,4:%i])\n",nFiles+3);
		fprintf(f,"print(pca)\n");
		fprintf(f,"plot(pca$x)\n");
		fclose(f);
	}

	writeLog("         ...OK\n");
	exit(0);
}

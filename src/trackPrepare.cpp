//============================================================================
// Name        : cpp.cpp
// Author      : a
// Version     :
// Copyright   : Your copyright notice
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "track_util.h"

//=====================================================================================
void bTrack::initProfile(){
	errStatus="init Profile";
	const char *memError="not enaugh memory. Try to increase  parameter \"step\"";

	getMem0(profile,profileLength, memError);					//== allocate memory
	if(profile==0) errorExit("%s",memError);
	if(strandFg){
		getMem0(profilec,profileLength, memError);			//== allocate memory for complement
		if(profilec==0) errorExit("%s",memError);
	}

	int zero=NA;
	for(int i=0; i<profileLength; i++){
		profile[i]=zero;								//== fill the profile with NA
		if(profilec) profilec[i]=zero;
	}
	errStatus=0;
};
//================================================================= Add segment
void bTrack::addSgm(ScoredRange *bed, float *prof){
	if(!checkRange(bed)) return;
	if(intervFlag) bed->score=1;
	else bed->score+=NA;
	int p1=pos2filePos(bed->chrom, bed->beg);
	int p2=pos2filePos(bed->chrom, bed->end);
	if(p1<0) p1=0; if(p2>=profileLength) p2=profileLength-1;
	if(p1==p2){
		prof[p1]+=bed->score*(bed->end-bed->beg);
	}
	else{
		int d=(int)(bed->score)*(stepSize*(p1 + 1 - curChrom->base) - bed->beg);
		prof[p1]+=d;
		for(int i=p1+1; i<p2; i++){
			prof[i]+=bed->score*stepSize;
		}
		d=(int)(bed->score)*(bed->end - (p2-curChrom->base)*stepSize);
		prof[p2]+=d;
	}
}

void bTrack::addSgm(char strnd, ScoredRange *bed, float *prof, float *profc, int strndFg){
	if(strndFg){								            //========== Strand dependent
		if(strnd=='-') addSgm(bed, profc);
		else 		   addSgm(bed, prof );
	}
	else			   addSgm(bed, prof);
}

//=====================================================================================
int getTypeByExt(const char *sext){
	char bext[80];
	strtoupper(strcpy(bext,sext));
	if(strcmp(bext,"BED")==0) return BED_TRACK;
	if(strcmp(bext,"WIG")==0) return WIG_TRACK;

	if(strcmp(bext,"BEDGRAPH"	)==0) return BED_GRAPH;
	if(strcmp(bext,"BED_GRAPH"	)==0) return BED_GRAPH;
	if(strcmp(bext,"B_GRAPH"	)==0) return BED_GRAPH;
	if(strcmp(bext,"BGR"		)==0) return BED_GRAPH;

	if(strcmp(bext,"B_PEAK"		)==0) return BROAD_PEAK;
	if(strcmp(bext,"BPEAK"		)==0) return BROAD_PEAK;
	if(strcmp(bext,"BROAD_PEAK"	)==0) return BROAD_PEAK;
	if(strcmp(bext,"BROADPEAK"	)==0) return BROAD_PEAK;

	if(strcmp(bext,"MODEL"	)==0) return MODEL_TRACK;
	if(strcmp(bext,"MOD"	)==0) return MODEL_TRACK;
	if(strcmp(bext,"MDL"	)==0) return MODEL_TRACK;

	return 0;
}
int getTrackType(const char *fname){
	const char *sext=getExt(fname);
	int tt=sext ? getTypeByExt(sext) : 0;
	return tt;
}
//============================================================ read track definition line
void bTrack::trackDef(char *s){
	char bb[100], *st;
	if(trackType==0){
		trackType=BED_TRACK;
		st=getAttr(s,"type",bb);
		if(st!=0){
			strtoupper(st);
			if(strncmp(st,"WIG",3)==0) 	    {trackType=WIG_TRACK; strandFg=0;}
			if(strncmp(st,"BEDGRAPH" ,8)==0) trackType=BED_GRAPH;
			if(strncmp(st,"BROADPEAK",8)==0) trackType=BROAD_PEAK;
		}
	}
	st=getAttr(s,"name", bb);
	if(st!=0) trackName=strdup(st);
}

//=====================================================================================
int countFields(char *b){
	int n=0;
	for(char *s=b; s; n++){
		char *ss=strchr(s,' ');
		if(ss==0) ss=strchr(s,'\t');
		if(ss) ss++;
		s=ss;
	}
	return n;
}

int checkWig(char *b){
	if(strncmp(b,(char *)"variableStep", 12)==0) return WIG_TRACK;
	if(strncmp(b,(char *)"fixedStep"   ,  9)==0) return WIG_TRACK;
	if(countFields(b)==4) return BED_GRAPH;
	return WIG_TRACK;
}


int ndeb=0;
void bTrack::readTrack(const char *fname, int cage){
	char *chrom=0;
    int fg=-1;
	int span=1, step=1;
	long start=0, beg=0, end=0;
	float score=0;
    ScoredRange bed;
    const int BUFFSIZE=20000;
	char buff[BUFFSIZE], abuf[1000], chBuf[100], *s;
	long i=0;
	char strand=0;
	int nStrand=0;
	strandFg=strandFg0;
	trackType=getTrackType(fname);
	intervFlag=(trackType == BED_TRACK) ? intervFlag0 : 0;
	lScale=intervFlag ? LIN_SCALE : logScale;

	FILE *f=xopen(fname, "rt"); setvbuf ( f , NULL , _IOFBF , 65536 );
	for(;(s=fgets(buff,BUFFSIZE,f))!=0; i++){
		strtok(s,"\n\r");
		int dataFg=1;									 	//======== remove end-of-line signes
		if(i%10000000 ==0)
			{verb("%li  %s...\n",i,curChrom->chrom); }
//if(i>10000000) return;
		if(EmptyString(buff)) continue;

		if(strncmp(s,"browser"  ,7)==0) continue;					//========= browser line => skip
		if(strncmp(s,"track"    ,5)==0) {trackDef(s); continue;}
		if(strncmp(s,"#bedGraph",9)==0) {trackType=BED_GRAPH; continue;}
		if(*s=='#' || *s==0) continue;							//======== comment line
		if(trackType==WIG_TRACK) trackType=checkWig(s);

		switch(trackType){
			case BED_TRACK:										//======== BED
				score=1;										//======== default score=1
				strand=0;										//======== default strand=unknown
				chrom=strtok(s," \t\n\r");						//======== find chrom field
				beg=atol(strtok(0," \t\r\n"));					//======== find beg field
				end=atol(strtok(0," \t\n\r"));					//======== find end field
				s=strtok(0," \t\n"); if(s==0) break;			//======== ignore name field
				s=strtok(0," \t\n"); if(s==0) break;			//======== take score field
				if(*s!=0 && (isdigit(*s) || *s=='-')) score=atof(s);
				if(intervFlag!=NONE) score=1;
				s=strtok(0," \t\n"); if(s==0) break;			//======== take strand field
				if(*s!=0 && *s!='.') {strand=*s; nStrand++;}
				if(intervFlag==GENE || intervFlag==NONE) break;
				if(intervFlag==GENE_BEG){
					if(strand=='-') beg=end-1;
					else 			end=beg+1;
					break;
					}
				if(intervFlag==GENE_END) {
					if(strand=='-') end=beg+1;
					else  			beg=end-1;
					break;
				}
				s=strtok(0,"\t"); if(s==0) break;			//======== skip tick1
				s=strtok(0,"\t"); if(s==0) break;			//======== skip tick2
				s=strtok(0,"\t"); if(s==0) break;			//======== skip rgb
				s=strtok(0,"\t"); if(s==0) break;			//======== number of exons
//======================= read gene structure
				{
				int nn=atoi(s);
				if(nn<2) break;
				int lExn[nn];
				int posExn[nn];
				bed.chrom=chrom; bed.score=score=1;
				for(int i=0; i<nn; i++) lExn[i]=atoi(s=strtok(0,",\t"));
				for(int i=0; i<nn; i++) posExn[i]=atoi(s=strtok(0,",\t"));
				int genePos=beg;
//===== now 'beg' and 'end' are positions relative to gene beg/
				if((intervFlag&EXON)==EXON){
					for(int i=0; i<nn; i++){
						end=(beg=posExn[i])+lExn[i];
						if(intervFlag==EXON_BEG){
							if(strand=='-') beg=end-1;
							else 			end=beg+1;
						}
						if(intervFlag==EXON_END){
							if(strand=='-') end=beg+1;
							else  			beg=end-1;
						}
					bed.beg=genePos+beg; bed.end=genePos+end;
						addSgm(strand, &bed, profile, profilec, strandFg);
					}
					break;
				}
				if((intervFlag&IVS)==IVS){
					for(int i=0; i<nn-1; i++){
						beg=posExn[i]+lExn[i];
						end=posExn[i+1];
						if(intervFlag==IVS_BEG){
							if(strand=='-') beg=end-1;
							else 			end=beg+1;
						}
						if(intervFlag==IVS_END) {
							if(strand=='-') end=beg+1;
							else  			beg=end-1;
						}
						bed.beg=genePos+beg; bed.end=genePos+end;
						addSgm(strand, &bed, profile, profilec, strandFg);
					}
					break;
				}
				}
				break;
			case BED_GRAPH:
				chrom=strtok(s," \t\n\r");						//======== find chrom field
				beg=atol(strtok(0," \t\r\n"));					//======== find beg field
				end=atol(strtok(0," \t\n\r"));					//======== find end field
				score=atof(strtok(0," \t\n\r"));				//======== find score field
				break;
			case WIG_TRACK:
				if(strncmp(s,"variableStep", 12)==0){			//======== read parameters for variableStep
					chrom=getAttr(buff,(char *)"chrom",chBuf);
					span=1; fg=0; dataFg=0;
					s=getAttr(buff,(char *)"span",abuf);
					if(s!=0) span=atoi(s);
				}
                else if(strncmp(s,(char *)"fixedStep", 9)==0){	//======== read parameters for fixedStep
                	fg=1; dataFg=0;
					chrom=getAttr(buff,(char *)"chrom",chBuf);
					s=getAttr(buff,(char *)"span",abuf);
					if(s!=0) span=atoi(s);
					s=getAttr(buff,(char *)"step",abuf);
					if(s!=0) step=atoi(s);
					s=getAttr(buff,(char *)"start",abuf);
					if(s!=0) start=atol(s);
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
					else errorExit("wrong WIG format");
				}
				break;
			case BROAD_PEAK:
				chrom=strtok(s," \t\n\r");						//======== find chrom field
				beg=atol(strtok(0," \t\r\n"));					//======== find beg field
				end=atol(strtok(0," \t\n\r"));					//======== find end field
				strtok(0," \t\n\r");							//======== skip name
				s=strtok(0," \t\n\r");
				if(bpType==BP_SCORE) score=atof(s);				//======== find score field
				s=strtok(0," \t\n"); if(s==0) break;			//======== take strand field
				if(*s!=0 && *s!='.') {strand=*s; nStrand++;}
				if(bpType==BP_SCORE) break;
				s=strtok(0," \t\n"); if(s==0) break;			//======== take signal field
				if(bpType==BP_SIGNAL) {score=atof(s); break;}
				s=strtok(0," \t\n"); if(s==0) break;			//======== take pval field
				if(bpType==BP_LOGPVAL) {score=atof(s); break;}
				break;
			default:
				errorExit("track type undefined or unknown"); break;
		}
		bed.chrom=chrom; bed.beg=beg; bed.end=end; bed.score=score;
		if(cage > 0) bed.end=bed.beg+cage;
		if(cage < 0) bed.beg=bed.end+cage;
		if(dataFg && intervFlag!=IVS && intervFlag!=EXON) {
			addSgm(strand, &bed, profile, profilec, strandFg);
		}
	}
	if(nStrand==0) strandFg=0;
	fclose(f);
}
//=====================================================================================
//========================================================================================
float bTrack::normProf(float x){
	float a=x/stepSize;
	if(lScale==LOG_SCALE)
		{if(a<=-1) a=-0.9999; a=a*scaleFactor; a=log(a+1);}
	return a;
}

void testDistrib(){
	float dd[1000]; memset(dd,0,sizeof(dd));
	for(int i=0; i<profileLength; i++){
		if(profile[i]==NA) dd[0]++;
		else{
			float x=log(1+profile[i]);
			int k=(int)(x/10*1000)+1;
//			if(k>=1000){
//				deb("overflow profile value  %f",x); k=999;
//			}
			dd[k]++;
		}
	}

	FILE *f=xopen("dstr","wt");
	for(int i=0; i<1000; i++){
		fprintf(f,"%5.2f\t%6f\n",i*10./1000.,dd[i]);
	}
	fclose(f);
}

void bTrack::finProfile(){
//	testDistrib();
	lScale=LOG_SCALE;
	scaleFactor=scaleFactor0;
	errStatus="finProfile";
	double lprof=0.;
	av=0.;            // Average profile value
	sd=1.;          // Standard deviation
	minP= 5.e+20;      // Minimal profile value
	maxP=-5.e+20;      // Maximal profile value
	//============================ calculate min, max, average, std deviation
	double x2=0;

	if(NAFlag==0){
		if(trackType==BED_TRACK){
			for(int i=0; i<profileLength; i++){
				if(profile [i] == NA) profile [i]=0;
				if(profilec && profilec[i] == NA) profilec[i]=0;
			}
		}
	}

	int nn=0;
	for(int i=0; i<profileLength; i++){
		if(profile[i] != NA){       // we take into account only valid profile values
			float z=profile[i]/stepSize;
			float a=normProf(profile[i]);
			av+=z; x2+=z*z; nn++;
			lprof+=1;
			if(a < minP) minP=a;
			if(a > maxP) maxP=a;
		}
		if(profilec && profilec[i] != NA){       // we take into account only valid profile values
			float z=profilec[i]/stepSize;
			float a=normProf(profilec[i]);
			av+=z; x2+=z*z; nn++;
			lprof+=1;
			if(a < minP) minP=a;
			if(a > maxP) maxP=a;
		}
	}
	av/=nn;
	if(maxP < minP){errorExit("=== !!!!  The profile contains no data:  min=%.2e max=%.2e !!!===\n",minP,maxP);}
	if(maxP==minP){	errorExit("=== !!!!  The profile contains only zero: min=%.2e max=%.2e  !!!===\n",minP,maxP);}
	x2-=av*av*nn;
	sd=sqrt(x2/(nn-1));
//	if(lScale==AUTO_SCALE  && maxP/sd > 10){
//		lScale=LOG_SCALE;
//		finProfile();
//		return;
//	}
//	else{
//		lScale=LIN_SCALE;
//	}
//
	for(int i=0; i<profileLength; i++){
		if(profile[i] != NA){       				// we take into account only valid profile values
			profile[i]=normProf(profile[i]);
		}
		if(profilec!=0 && profilec[i] != NA){       // we take into account only valid profile values
			profilec[i]=normProf(profilec[i]);
		}
	}

	float bScale=250./(maxP-minP);           		 //======================== scale coefficient
	getMem0(byteprofile,profileLength, "bTrack::finProfile #1");
	getMem0(byteprofilec,profileLength, "bTrack::finProfile #2");

	for(int i=0; i<profileLength; i++){
		if(profile[i] !=NA)  {
			byteprofile[i] =int(1.+bScale*(profile[i]-minP));
		}
		else                 byteprofile[i] =0;       // undefined values are presented as 0 in the byte profile
		if(profilec!=0 && profilec[i] !=NA){
			byteprofilec[i]=int(1.+bScale*(profilec[i]-minP));
		}
		else                 byteprofilec[i]=0;       // undefined values are presented as 0 in the byte profile
	};
	errStatus=0;
};
//===============================================================================================
//=================================================================================
//#        Write the byte profile in two files.
//#            *.prm - file with average min and max values
//#            *.bprof - byte profile

//=================================================================================
unsigned int getChkSum(unsigned char *b, int n){
	unsigned int chksum=0;
	for(int i=0; i<n; i++){
		chksum+=b[i];
	}
	return chksum;
}

void bTrack::writeProfilePrm(){
    //============================================= Write parameters
	char prmFname[4096];
	makeFileName(prmFname,profPath,name,PRM_EXT);
	verb("write prm %s\n",prmFname);
    FILE *f=xopen(prmFname, "wt"); if(f==0)return;


	fprintf(f,"#====== THIS IS GENERATED FILE. DO NOT EDIT!  ========\n");
	fprintf(f,"version=%s\n",version);
	fprintf(f,"#=== Input data ===\n");

    fprintf(f,"trackType=%i\n",trackType);
	if(trackType == MODEL_TRACK)
		fprintf(f, "input=%s\n", model.definition);
	else
		fprintf(f,"input=%s\n",name);

    fprintf(f,"#=== Parameters ===\n");
    fprintf(f,"step=%i\n",stepSize);
    fprintf(f,"strand=%i\n",strandFg);
    if(trackType == BROAD_PEAK)  fprintf(f,"bpType=%i\n",bpType);
    if(trackType == BED_TRACK)   fprintf(f, "intervFlag=%i\n", intervFlag);
    fprintf(f,"NA=%i\n",NAFlag);
    if(cage) fprintf(f, "cage=%i\n", cage);
    fprintf(f,"scaleFactor=%f\n",scaleFactor);

    fprintf(f,"#=== Statistics ===\n");
    fprintf(f,"min=%g\n",minP);
    fprintf(f,"max=%g\n",maxP);
    fprintf(f,"average=%g\n",av);
    fprintf(f,"stdDev=%g\n",sd);
    unsigned int chk=getChkSum(byteprofile, profileLength);
    if(byteprofilec && strandFg) chk+=getChkSum(byteprofilec, profileLength);
    fprintf(f,"chksum=%08x\n",chk);

    fprintf(f,"#=== Scale ===\n");
    fprintf(f,"lscale=%i\n",lScale);
    //==== calculate distribution
    int dstr[256];
    float cdstr[256],nn=0;
    memset(dstr,0,sizeof(dstr));
    for(int i=0; i<profileLength; i++){
    	dstr[byteprofile[i]]++;
    	if(byteprofilec && strandFg) dstr[byteprofilec[i]]++;
    }

    cdstr[0]=0; nn=0;
    for(int i=0; i<256; i++){
    	nn+=dstr[i];
    	if(i>0) cdstr[i]=nn;
    }
    fprintf(f,"\n#***** Distribution ***** nn=%i  pLength=%i\n",(int)nn,profileLength);
    for(int i=0; i<256; i++){
    	fprintf(f,"#%3i\t%i\t%.0f\t%5.2f\t%5.2f\n",i,dstr[i],cdstr[i],
    	                                            cdstr[i]*100./nn, (nn-cdstr[i])*100./nn);
    }

	fclose(f);
}


void bTrack::writeByteProfile(){
	//============================================= Write bytes
	char binFname[4096]; makeFileName(binFname,profPath,name,BPROF_EXT);
	verb("write bprof %s\n",binFname);
	FILE* f=xopen(binFname, "wb"); if(f==0)return;
	fwrite(byteprofile,profileLength,1,f);
	if(byteprofilec && strandFg)
		fwrite(byteprofilec,profileLength,1,f);
	fclose(f);
	return;
}
//=================================================================================


void bTrack::makeBinTrack(){
	verb ("******   Make binary track <%s>   ******\n", name);
	writeLog("    Make binary track <%s> ",name);
	//===================================================== prepare track file name
	if (pcorProfile!=0)	 verb("===          pcorProfile=  <%s>\n",pcorProfile);
	trackType=getTrackType(name);
	//=====================================================================================
	intervFlag=(trackType==BED_TRACK) ? intervFlag0 : 0;

	if(trackType==MODEL_TRACK) {
		model.readMap(name);
		model.create();
		return;
	}

	initProfile();                   //============ Allocate arrays
	char pfil[4096];
	readTrack(makeFileName(pfil,trackPath,name));		//============ Read tracks
	//======================================================================
	verb("Finalize profiles... \n");
	finProfile(); //============ Calculate min,max,average; convert to bytes
	verb("Write profiles... \n");
	writeProfilePrm();
	writeByteProfile();					//============ write binary profiles
//	if(lAuto) trackAutoCorrelation();
	writeLog("         ...OK\n");

	return;
}


//================================================================================
Model::Model(){
	definition=0; nTerm=0; trackType=MODEL_TRACK;
}

//================================================================================
void Model::readMap(char *fnam){
	name=strdup(fnam);
	char modName[4096]; makeFileName(modName, trackPath, name);
	FILE * f=xopen(modName,"rt");
	char b[4096], *s;
	for(;(s=fgets(b,sizeof(b),f))!=0;){
		strtok(b,"\r\n#");
		s=skipSpace(b);
		if(*s=='#') continue;
		if(strlen(s)==0) continue;
		if(trm[nTerm].read(s)) nTerm++;
	}
	fclose(f);
	char bb[1024];*b=0;
	for(int i=0; i<nTerm; i++){
		sprintf(bb,"%s(%f)*[%s]",i ? "+" : "" ,trm[i].mult,trm[i].fname);
		strcat(b,bb);
	}
	definition=strdup(b);
}
//================================================================================
void delTailSpase(char *s){
	for(char *sx=s+strlen(s)-1; isspace(*sx); sx--) *sx=0;
}
int Term::read(char *s){
	char *s1=strtok(s,"*"), *s2=strtok(0,"*");
	delTailSpase(s1);
	if(s2==0 || *s2==0) {
		mult=1; fname=strdup(s1);
	}
	else{
		mult=atof(s1);
		delTailSpase(s2); s2=skipSpace(s2);
		fname=strdup(s2);
	}
	return 1;
}

//===================================================================================
bTrack tmpBTrack;
float* tmpProf;

void Term::make(){
	tmpBTrack.name=fname;
	if(!tmpBTrack.check(fname))	tmpBTrack.makeBinTrack();
}

void Term::add(){
//	char b[4096];
	tmpBTrack.read(fname);
	verb("make Model. Add <%s>\n",fname);
	for(int i=0; i<profileLength; i++){
		float x=tmpBTrack.getVal(tmpBTrack.bytes[i]);
		if(tmpBTrack.hasCompl) x+=tmpBTrack.getVal(tmpBTrack.cbytes[i]);
		x*=mult;
		profile[i]=x;
	}
	tmpBTrack.clear();
}
//=====================================================================================
void Model::create(){
	int scale=logScale; logScale=LIN_SCALE;
	strandFg=0;

	for(int i=0; i<nTerm; i++)	trm[i].make();

	initProfile();
	for(int i=0; i<nTerm; i++)	trm[i].add();

	finProfile();
	writeProfilePrm();
	writeByteProfile();
	logScale=scale;
}

//=======================================================================================
//=======================================================================================
void CageMin(const char* fname1, const char * fname2){	//read Cage files
	float *pr1;
	char b[4096], b1[4096], b2[4096], *s;
	strcpy(b1,fname1); s=strrchr(b1,'.'); if(s) *s=0;
	strcpy(b2,fname2); s=strrchr(b2,'.'); if(s) *s=0;
	sprintf(b,"MIN_%s_%s.",b1,b2);

	bTrack1.name=strdup(b);
	bTrack1.initProfile();
	fname1=makeFileName(b,trackPath,fname1);
	bTrack1.readTrack(fname1,cage);
	bTrack1.trackType=CAGE_MIN;
	pr1=profile; profile=0;;
	bTrack2.initProfile();
	fname2=makeFileName(b,trackPath,fname2);
	bTrack2.readTrack(fname2,-cage);bTrack2.name=strdup(fname2);

	verb("read OK\n");
	for(int i=0; i<profileLength; i++){
		profile[i]=(profile[i]<pr1[i]) ? profile[i] : pr1[i];
	}
	verb("min OK\n");
	bTrack1.finProfile();
	verb("fin OK\n");
	bTrack1.writeProfilePrm();
	verb("prm OK\n");
	bTrack1.writeByteProfile();
	verb("!!! OK\n");
}




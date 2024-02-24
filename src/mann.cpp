/*
 * mann.cpp
 * Mann-whitney test
 *  Created on: Feb 26, 2013
 *      Author: Mironov
 */


#include "track_util.h"


//==================================================== Testing
const int nSet1=200;
const int nSet2=100;
double set1[nSet1], set2[nSet2];
void makeSet(int nSet, double*set, double e, double sigma){
	for(int i=0; i<nSet; i++)
		set[i]=rGauss(e,sigma);
}
void makeSets(){
	makeSet(nSet1, set1, 0.,1);
	makeSet(nSet2, set2, 0.1,1);
}
void writeSet(double *set, int nSet, const char*fname){
	FILE *f=xopen(fname,"wt");
	for(int i=0; i<nSet; i++)
		fprintf(f,"%f\n",set[i]);
	fclose(f);
}


//==================================================== Normal distribution
int stdLength=71;
double stdNormal[71] = {
                       5.000000E-01, 4.601721E-01, 4.207403E-01, 3.820886E-01, 3.445783E-01, //0.0-0.4
                       3.085375E-01, 2.742531E-01, 2.419636E-01, 2.118553E-01, 1.840601E-01, //0.4-0.9
                       1.586553E-01, 1.356661E-01, 1.150697E-01, 9.680055E-02, 8.075671E-02, //1.0-1.4
                       6.680723E-02, 5.479929E-02, 4.456543E-02, 3.593027E-02, 2.871649E-02, //1.4-1.9
                       2.275006E-02, 1.786436E-02, 1.390340E-02, 1.072408E-02, 8.197529E-03, //2.0-2.4
                       6.209680E-03, 4.661222E-03, 3.467023E-03, 2.555191E-03, 1.865880E-03, //2.4-2.9
                       1.349967E-03, 9.676712E-04, 6.872021E-04, 4.834825E-04, 3.369808E-04, //3.0-3.4
                       2.326734E-04, 1.591457E-04, 1.078301E-04, 7.237243E-05, 4.811552E-05, //3.4-3.9
                       3.168603E-05, 2.066872E-05, 1.335410E-05, 8.546021E-06, 5.416953E-06, //4.0-4.4
                       3.400803E-06, 2.114643E-06, 1.302316E-06, 7.943527E-07, 4.798695E-07, //4.4-4.9
                       2.871050E-07, 1.701223E-07, 9.983440E-08, 5.802207E-08, 3.339612E-08, //5.0-5.4
                       1.903640E-08, 1.074622E-08, 6.007653E-09, 3.326052E-09, 1.823579E-09, //5.4-5.9
                       9.901219E-10, 5.323753E-10, 2.834714E-10, 1.494721E-10, 7.804901E-11, //6.0-6.4
                       4.035794E-11, 2.066525E-11, 1.047862E-11, 5.261569E-12, 2.616130E-12, //6.4-6.9
                       1.288081E-12, //7.0
 };


double pisq = sqrt(2 * PI);
//==================================================== Normal density
double normalDensity(double x) {
    return exp( -x * x / 2) / pisq;
  }
//=================================================== Nrmal distribution ( 1-F(x) )
double normalProbGrater(double x) {
    if (x == 0)return 0.5;
    if (x < 0)return 1 - normalProbGrater(-x);
    int idx = (int) ( (x + 0.00001) / 0.1);
    if (idx < 0)return 0.5;
    double w = 0;
    if (idx >= stdLength - 1) {
      double z = (1 - 1 / (x * x));
      w = normalDensity(x) / x * z;				//== assimptotic
    }
    else {
      double k = (x - idx * 0.1) / 0.1;			//== interpolation
      w = stdNormal[idx] * (1 - k) + stdNormal[idx + 1] * k;
    }
    if (w <= 0) w = 0;
    return w;
  }




//================================================ comparator for sort double
int doubleCmp(const void *dd1, const void *dd2){
	double *d1=(double *)dd1;
	double *d2=(double *)dd2;
	if(*d1>*d2) return  1;
	if(*d1<*d2) return -1;
	return 0;
}


//double *sample(double *set, int nset,int nSample){
//	double *subset;
//	errStatus="mann subset";
//	getMem(subset,nSample);
//	for(int i=0; i<nSample; i++){
//		int ii=randInt(nset);
//		subset[i]=set[ii];
//	}
//	errStatus=0;
//	return subset;
//}
statTest *MannWhitney0( double *set1, int nSet1, double *set2,int nSet2);
//================================================ Mann-Whitney test
statTest *MannWhitney( double *set1, int nSet1, double *set2,int nSet2){
	return MannWhitney0(set1,nSet1,set2,nSet2);
}
statTest *MannWhitney0( double *set1, int nSet1, double *set2,int nSet2){
	if(nSet1==0) {writeLog("No background data"); return 0;}
	if(nSet2==0) {writeLog("No foreground data"); return 0;}
	qsort((void *)set1,nSet1,sizeof(double),doubleCmp);
	qsort((void *)set2,nSet2,sizeof(double),doubleCmp);


	double u=0;								//=== Mann-Whithey statistics
	for(int i1=0, i2=0; i1<nSet1; i1++){
		for(;set1[i1] > set2[i2] && i2<nSet2; i2++);
		u+=i2;
	}
	double dSet1=nSet1, dSet2=nSet2;		//=== Calculation should be done in double
	double e=0.5*(dSet1*dSet2);							//=== Mean
	double sigma=sqrt(dSet1*dSet2*(dSet1+dSet2+1)/12);	//=== standatd deviation
	double z=(u-e)/sigma;								//=== Z-score
	double zz=(z > 0)? z:-z;							//=== positive Z-score
	double p_value=normalProbGrater(zz)*2;				//=== normal approximation
	statTest *st=new statTest();						//=== store statistical test
	st->n1=nSet1;
	st->n2=nSet2;
	st->u =u;
	st->e=e;
	st->sigma=sigma;
	st->z=z;
	st->pVal=p_value;
	return st;
}




//======================================================= main for testing
int zxmain(int argc, char *argv[]){
	makeSets();
	writeSet(set1,nSet1,"set1");
	writeSet(set2,nSet2,"set2");
	MannWhitney(set1, nSet1, set2, nSet2);
	return 0;
}



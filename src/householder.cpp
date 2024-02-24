//https://rosettacode.org/wiki/QR_decomposition#C


#include "track_util.h"


//=========================================
Matrix::Matrix(int nn, double *a){
	init(nn,a);
}
//=========================================
Matrix::Matrix(int mm){
	init(mm);
}
//=========================================
Matrix::Matrix(Matrix *mtx){
	init(mtx->n,mtx->values);
}
//=========================================
void Matrix::init(int nn){
	n=nn;
	getMem(values,n*n,"aaaa"); zeroMem(values,n*n);
}


//=========================================
void Matrix::init(int nn, double *a){
	init(nn); set(a);
}


//=========================================
void Matrix::transpose(){
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			double t = get(i,j);
			set(i,j,get(j,i));
			set(j,i,t);
		}
	}
}
//==========================================
void Matrix::printMtx(FILE *f){
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(f," %8.3f", get(i,j));
		}
		fprintf(f,"\n");
	}
	fprintf(f,"\n");
}
//==========================================
void Matrix::printMtx(){
	printMtx(stdout);
}
//=========================================


void mult(Matrix *res, Matrix *x, Matrix *y){
	double z=0;
	for (int i = 0; i < x->n; i++){
		for (int j = 0; j < y->n; j++){
			z=0;
			for (int k = 0; k < x->n; k++){
				z+= x->get(i,k) * y->get(k,j);
			}
			res->set(i,j,z);
		}
	}
}
//=================  Multiply matrices ============
Matrix *mult(Matrix *x, Matrix *y){
	if(x->n!=y->n) return 0;
	Matrix *res=new Matrix(x->n);
	mult(res,x,y);
	return res;
}


//=================  matrix minor ============
Matrix *minorMtx(Matrix *mtx, int d){
	Matrix* m = new Matrix(mtx->n);
	for (int i = 0; i < d; i++)
			m->set(i,i,1);
	for (int i = d; i < mtx->n; i++)
		for (int j = d; j < mtx->n; j++)
			m->set(i,j,mtx->get(i,j));
	return m;
}


//=========== take c-th column of m, put in v
double* mcol(Matrix *mtx, double *v, int c){
	for (int i = 0; i < mtx->n; i++)
		v[i] = mtx->get(i,c);
	return v;
}


//=============================== m = I - v v^T
Matrix *vmul(double v[], int n){
	Matrix *xxx = new Matrix(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			xxx->set(i,j,-2 *  v[i] * v[j]);
	for (int i = 0; i < n; i++)
		xxx->set(i,i,xxx->get(i,i)+1);
	return xxx;
}


//============================== c = a + b * s
double *vecMultAdd(double a[], double b[], double s, double c[], int n){
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
}


//=============================== ||x||
double vecNorm(double x[], int n){
	double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}
 
//================================ y = x / d
double* vecDiv(double x[], double d, double y[], int n){
	if(d==0) return y;
	for (int i = 0; i < n; i++) y[i] = x[i] / d;
	return y;
}


//========================================================
void householder(Matrix *mtx, Matrix **R, Matrix **Q){
	int mm=mtx->n;
	Matrix *q[mm];
	Matrix *z = mtx, *z1;
	for (int k = 0; k < mtx->n && k < mtx->n - 1; k++) {
		double e[mm], x[mm], a;
		z1 = minorMtx(z, k);
		if (z != mtx) {delete z; z=0;}
		z = z1;
		mcol(z, x, k);
		a = vecNorm(x, mm);
		if (mtx->get(k,k) > 0) a = -a;
 
		for (int i = 0; i < mm; i++)
			e[i] = (i == k) ? 1 : 0;
		vecMultAdd(x, e, a, e, mm);
		vecDiv(e, vecNorm(e, mm), e, mm);


		q[k] = vmul(e, mm);
		z1 = mult(q[k], z);
		if (z != mtx) {delete z; z=0;}
		z = z1;
	}
	if (z != mtx) {delete z; z=0;}
	*Q = q[0];
	*R = mult (q[0], mtx);
	for (int i = 1; i < mtx->n && i < mtx->n - 1; i++) {
		z1 = mult(q[i], *Q);
		if (i > 1) {delete *Q; Q=0;}
		*Q = z1;
		delete q[i]; q[i]=0;
	}
	delete q[0]; q[0]=0;
	z = mult(*Q, mtx);
	delete *R; *R=0;
	*R = z;
	(*Q)->transpose();
}
//=========================================================
//===============  Eigen Vectors ==========================
Matrix * eigenVectors(Matrix *x, double *EValues, int nIter, double precsision){
	Matrix *R=0;
	Matrix *Q=0;
	Matrix *EVal=new Matrix(x);
	Matrix *Vect= minorMtx(EVal,EVal->n);
	double di=0, noDi=0;
	int iter=0;
	for(; iter<nIter; iter++){
		householder(EVal, &R, &Q);
		mult(EVal,R,Q);
		Matrix *xqz=Vect;
		Vect=mult(Vect,Q);
		delete xqz; xqz=0;
		for(int i=0; i<x->n; i++){
			for(int j=0; j<x->n; j++){
				double w=EVal->get(i,j);
				w*=w;
				if(i==j) di+=w; else noDi+=w;
			}
		}
		if(noDi/di < precsision) break;
	}
	if(EValues){
		for(int i=0; i<EVal->n; i++)
			EValues[i]=EVal->get(i,i);
	}
	delete R; R=0;
	delete Q; Q=0;
	delete EVal; EVal=0;


return Vect;
}


//=========================================================
//===============  Eigen Vectors ==========================
double in[] = {
	 6   ,4  ,1.2 , 1,
	 4   ,7  ,2   , 1,
	 1.2 ,2  ,8   ,-1,
	 1   ,1  ,-1  , 9,
	};


void hh_test(){
	int n=4;
	int nIter=200;
//	for(int i=0; i<n*n; i++) in[i]/=10;
	Matrix *x = new Matrix(n,in);


	puts("input"); x->printMtx();


	double eValues[n];
	Matrix *Vect=eigenVectors(x,eValues,nIter,1.e-3);


	puts("Vect"); Vect->printMtx();
	puts("EVal");
	for(int i=0; i<n; i++)
		printf("%.5f; ",eValues[i]);
	printf("\n");
 
	del(x);
	exit(0);
	}


/*
Output:
Q
    0.846   -0.391    0.343    0.082    0.078
    0.423    0.904   -0.029    0.026    0.045
   -0.282    0.170    0.933   -0.047   -0.137
   -0.071    0.014   -0.001    0.980   -0.184
    0.141   -0.017   -0.106   -0.171   -0.969


R
   14.177   20.667  -13.402
   -0.000  175.043  -70.080
    0.000    0.000  -35.202
   -0.000   -0.000   -0.000
    0.000    0.000   -0.000


Q * R
   12.000  -51.000    4.000
    6.000  167.000  -68.000
   -4.000   24.000  -41.000
   -1.000    1.000   -0.000
    2.000   -0.000    3.000
    */

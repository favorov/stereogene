/*
 * kernel.cpp
 *
 *  Created on: Aug 30, 2013
 *      Author: mironov
 */
#include "track_util.h"

Fourier::Fourier(int n){init(n);}
Fourier::Fourier(){re=im=datRe=datIm=0; length=0; err=0; re0=im0=0;}
double Complex::Mod(){return sqrt(re*re+im*im);}
Complex Complex::scalar(Complex otherC){
	Complex res=Complex();
	res.re = re * otherC.re + im * otherC.im;
	res.im = - re*otherC.im + im*otherC.re;
	return res;
}

void Fourier::init(int len){
	length=len; err=0;
	errStatus="Fourier init";
	getMem(re,len, "Fourier init #1");
	getMem(im,len, "Fourier init #2");
	getMem(datIm,len, "Fourier init #3");
	zeroMem(datIm,len);
	errStatus=0;
}

Fourier::~Fourier(){
	freeMem();
}

void Fourier::freeMem(){
	if(re) free(re); re=0;
	if(im) free(im); im=0;
	if(datIm) free(datIm); datIm=0;
}

void Fourier::setDat(double *reD){
	datRe=reD;
	zeroMem(datIm,length);
}
void Fourier::setDat(double *reD, double *imD){
	datRe=reD;
	datIm=imD;
}
void Fourier::calc(double *reD, double *imD, int deriv){
	double *dim0=datIm;
	setDat(reD,imD);
	calc(deriv);
	datIm=dim0;
}
void Fourier::calc(double *dRe, int deriv){
	setDat(dRe);
	calc(deriv);
}
void Fourier::calc0(double *dRe, int deriv){
	calc(dRe, deriv);
	re0=re[0]; im0=im[0];
	re[0]=0; im[0]=0;
}
void Fourier::calc(int deriv){
	if(err) return;
	fftl(length,datRe, datIm,re,im);


	for(int i=0; i<deriv; i++) derivat();
}
void Fourier::derivat(){
	double piL=PI*2/length;
	for(int i=0; i<length; i++){
		double rre=-i*piL*im[i], iim=i*piL*re[i];
		re[i]=rre; im[i]=iim;
	}
}

int corrFunc(double *x, double *y, double *rc, int l){
	if(norm(x,l)==0) return 0;
	if(norm(y,l)==0) return 0;
	Fourier frx(l); frx.calc(x,0);
	Fourier fry(l); fry.calc(y,0);
	double *ic; getMem(ic,l, "corrFunc #1");
	for(int i=0; i<l; i++){
		x[i]=frx.re[i]*frx.re[i]+frx.im[i]*frx.im[i];
		y[i]=fry.re[i]*fry.re[i]+fry.im[i]*fry.im[i];
		rc[i]=fry.re[i]*frx.re[i]+fry.im[i]*frx.im[i];
		ic[i]=frx.re[i]*fry.im[i]-fry.re[i]*frx.im[i];
	}
	frx.calc(x,0);
	fry.calc(y,0);
	for(int i=0; i<l; i++) {x[i]=frx.re[i]/l/l; y[i]=fry.re[i]/l/l;}
	frx.calc(rc,ic,0);
	for(int i=0; i<l; i++) {rc[i]=frx.re[i]/l/l;}
	return 1;
}


void Kernel::init(int n){
	length=n;
	errStatus="kernel init";
	getMem(kern,length, "kernel init #1");
	getMem(ckern,length, "kernel init #2");
	ft.init(n);
	cft.init(n);
	fx.init(n);
	fy.init(n);
}

void Kernel::fftx(double* x, int deriv){
	fx.calc0(x,deriv);
}
void Kernel::ffty(double* y, int deriv){
	fy.calc0(y,deriv);
}


void Kernel::fft(){
	ft.calc0(kern,0);
	cft.calc0(ckern,0);
}

double Kernel::scalar(Fourier *f1, Fourier *f2, Complex *c, bool complem){
	if(f1->length !=length) return 0;
	if(f2->length !=length) return 0;
	double re=0, im=0;
	Fourier *zft=complem ? &cft : &ft;
	for(int i=1; i<length/2; i++){
		double RaRb_plus_IaIb =f1->re[i]*f2->re[i] + f1->im[i]*f2->im[i]; //==== Re(a)*Re(b)+Im(a)*Im(b)
		double RaIb_minus_IaRb=f1->re[i]*f2->im[i] - f1->im[i]*f2->re[i]; //==== Re(a)*Im(b)-Im(a)*Re(b)
		//==  Re=SUM Re(kern)(Re(a)*Re(b)+Im(a)*Im(b)) + Im(kern)*Re(a)*Im(b)-Im(a)*Re(b)
		re+=zft->re[i]*RaRb_plus_IaIb + zft->im[i]*RaIb_minus_IaRb;
		//==  Im=SUM Im(kern)(Re(a)*Re(b)+Im(a)*Im(b)) - Re(kern)*Re(a)*Im(b)-Im(a)*Re(b)
		im+=zft->im[i]*RaRb_plus_IaIb - zft->re[i]*RaIb_minus_IaRb;
	}
	if(c!=0) {c->re=re; c->im=im;}
	return re;
}

double Kernel::dist(bool complem){
		return dist(&fx,&fy, complem);
}


double Kernel::dist(Fourier *f1, Fourier *f2, bool complem){
	Complex c0=Complex(), c1=Complex(), c2=Complex();

	double d12=scalar(f1,f2, &c0, complem);
	double d11=scalar(f1,f1, &c1, complem);
	double d22=scalar(f2,f2, &c2, complem);
	d11=c1.Mod(); d22=c2.Mod();
	if(d11==0 || d22==0) return -400;
	Fourier *zft=complem ? &cft : &ft;
	double dd;

	dd=f1->re0*f1->re0*zft->re0; prod11+=d11+dd; sprod11+=dd;
	dd=f1->re0*f2->re0*zft->re0; prod12+=d12+dd; sprod12+=dd;
	dd=f2->re0*f2->re0*zft->re0; prod22+=d22+dd; sprod22+=dd;

	nProd++;
	double cc=d12/sqrt(d11*d22);

	return cc;
}

double Kernel::dist(Fourier *f1, Fourier *f2, Fourier *fpc, bool complem){
	Complex c0=Complex(), c1=Complex(), c2=Complex();

	double d12=scalar(f1,f2, &c0, complem);//<xy>
	double d11=scalar(f1,f1, &c1, complem);//<xx>
	double d22=scalar(f2,f2, &c2, complem);//<yy>

	d11=c1.Mod(); d22=c2.Mod(); //|x|, |y|

	//partial correlation
	Complex cx=Complex(), cy=Complex(), cz=Complex();
	double dx = scalar(f1,fpc, &cx, complem);  //<xz>
	double dy = scalar(f2,fpc, &cy, complem);  //<xz>
	double dz2 = scalar(fpc,fpc, &cz, complem);//<zz>

	dz2=cz.Mod();
	dx = cx.Mod();
	dy = cy.Mod();

	double dxx = d11 - dx*dx/dz2;
	double dyy = d22 - dy*dy/dz2;
	double dxy = (cx.scalar(cy)).re; //<xz><yz>

	if(d11==0 || d22==0 || dz2==0) return -400;
	Fourier *zft=complem ? &cft : &ft;

	double dd;

	dd=f1->re0*f1->re0*zft->re0; prod11+=d11+dd; sprod11+=dd;
	dd=f1->re0*f2->re0*zft->re0; prod12+=d12+dd; sprod12+=dd;
	dd=f2->re0*f2->re0*zft->re0; prod22+=d22+dd; sprod22+=dd;
	nProd++;
//-------------------------------------
	return (d12 - dxy/dz2) /sqrt(dxx*dyy);	// (<xy>-<xz><yz>/z^2)/sqrt((x^2-|xz|/z^2)(y^2-|yz|/z^2))
}

void Kernel::makeKernel(int n){
	init(n);
	double d=0;
	for(int i=0; i<length; i++){
		int x1=i, x2=i-length;
		ckern[length-i-1]=kern[i]=kernVal(x1)+kernVal(x2);
		d+=kern[i];
	}
	for(int i=0; i<length; i++){
		ckern[length-i-1]/=d;
		kern[i]/=d;
	}
	fft();
}


NormKernel::NormKernel(){sigma=1; e=0; name= strdup("Normal_Kernel");}
NormKernel::NormKernel(double ee,double sgm, int l){
	sigma=sgm; e=ee; hasCompl=(e!=0);
	makeKernel(l);name= strdup("Normal_Kernel");
}

double Kernel::NSCorrection(double x, double val){
	double x0=x*stepSize;
	if(x0 < 100 && x0>-100){
		val*=(1-kernelNS);
	}
	return val;
}

double NormKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma;
	double val=exp(-x*x);
	return NSCorrection(x0,val);
}

LeftExpKernel::LeftExpKernel(){sigma=1; e=0; name= strdup("Left_Exp_Kernel");}
LeftExpKernel::LeftExpKernel(double ee,double sgm, int l){
	sigma=sgm; e=ee;
	makeKernel(l); name= strdup("Left_Exp_Kernel");
}
double LeftExpKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma; hasCompl=true;
	return NSCorrection(x0,(x>0)?0 : exp(x));
}


RightExpKernel::RightExpKernel(){sigma=1;e=0;name= strdup("Right_Exp_Kernel");}
RightExpKernel::RightExpKernel(double ee,double sgm, int l){
	sigma=sgm; e=ee;
	makeKernel(l);name= strdup("Right_Exp_Kernel");
}
double RightExpKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma; hasCompl=true;
	return NSCorrection(x0,(x<0)?0:exp(-x));
}




/*
 * kernel.cpp
 *
 *  Created on: Aug 30, 2013
 *      Author: Mironov
 */
#include "track_util.h"


double Complex::Mod(){return sqrt(re*re+im*im);}


Complex Complex::scalar(Complex otherC){
	Complex res=Complex();
	res.re = re * otherC.re + im * otherC.im;
	res.im = - re*otherC.im + im*otherC.re;
	return res;
}


Fourier::Fourier(int n){init(n);}


Fourier::Fourier(){re=im=datRe=datIm=0; length=0; err=0; spectrum=0; autocorr=0;}


void Fourier::init(int len){
	if(length==len) return;
	length=len; err=0;
	errStatus="Fourier init";
	if(re   ) xfree(re   , "free Fourier #1");
	if(im   ) xfree(im   , "free Fourier #2");
	if(datRe) xfree(datRe, "free Fourier #3");
	if(datIm) xfree(datIm, "free Fourier #4");
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
	if(re!=0) xfree(re,"re");
	if(im!=0) xfree(im,"im");
	if(datIm) xfree(datIm,"datIm");
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
	double *dIm0=datIm;
	setDat(reD,imD);
	calc(deriv);
	datIm=dIm0;
}
void Fourier::calc(double *dRe, int deriv){
	setDat(dRe);
	calc0(deriv);
}


void Fourier::norm(){
	for(int i=0; i<length; i++) {re[i]/=length; im[i]/=length;}
}


void Fourier::calc(int deriv){
	if(err) return;
	zeroMem(re,length);	zeroMem(im,length);
	calc(datRe, datIm,re,im);
	derivat(deriv);
}
void Fourier::calc0(int deriv){
	if(err) return;
	zeroMem(datIm,length);
	calc(deriv);
}




void Fourier::calc(double *dRe,double *dIm,double *rRe,double *rIm){
	fftl(length, dRe, dIm, rRe, rIm);
}


void Fourier::derivat(int deriv){
	for(int i=0; i<deriv; i++) derivat();
}




void Fourier::derivat(){
	double piL=PI*2/length;
	for(int i=0; i<length; i++){
		double rre=-i*piL*im[i], iim=i*piL*re[i];
		re[i]=rre; im[i]=iim;
	}
}


float *Fourier::getSpectrum(){
	getMem(spectrum,length,"Spectrum");
	for(int i=0; i<length; i++)
		spectrum[i]=(float)(re[i]*re[i]+im[i]*im[i]);
	return spectrum;
}


//================== calculate autocorrelations ====================
double *tmpDRe=0, *tmpDIm=0, *tmpIm=0;


double *Fourier::getAutoCorr(){
	getMem0(autocorr,length,"AutoCorr #1");
	getMem0(tmpDRe,length,"AutoCorr #2");
	getMem0(tmpDIm,length,"AutoCorr #3");
	getMem0(tmpIm,length,"AutoCorr #4");


	tmpDRe[0]=0;
	for(int i=1; i<length; i++){
		tmpDRe[i]=re[i]*re[i]+im[i]*im[i];
	}
	zeroMem(tmpDIm,length);
	calc(tmpDRe, tmpDIm, autocorr, tmpIm);
	return autocorr;
}


//=================================================================
//============================   Kernel ===========================
//=================================================================
Kernel *MakeKernel(int l){
	Kernel *kern=0;
	switch(kernelType){
	case KERN_NORM:
		kern=new NormKernel    (kernelProfShift, kernelProfSigma, l); break;
	case KERN_LEFT_EXP:
		kern=new LeftExpKernel (kernelProfShift, kernelProfSigma, l); break;
	case KERN_RIGHT_EXP:
		kern=new RightExpKernel(kernelProfShift, kernelProfSigma, l); break;
	case KERN_CUSTOM:
		kern=new CustKernel(kernelProfShift, kernelProfSigma, l); break;
	default: errorExit("Kernel not defined");  break;
	}
	return kern;
}


void Kernel::init(int n){
	length=n;
	errStatus="kernel init";
	if(kern ) xfree(kern ,"free Kern #1");
	if(ckern) xfree(ckern,"free Kern #1");
	kern=0; ckern=0;
	getMem(kern,length, "kernel init #1");
	getMem(ckern,length, "kernel init #2");
	ft.init(n);
	cft.init(n);
	fx.init(n);
	fy.init(n);
}


//============================== Calculate FFT
void Kernel::fftx(double* x, int deriv){
	fx.calc(x,deriv);
}
void Kernel::ffty(double* y, int deriv){fy.calc(y,deriv);}




//============================== Calculate FFT for direct & compl
void Kernel::fft(){
	ft.calc(kern,0);
	cft.calc(ckern,0);
}


//============================== Kerneled scalar prod =\int f(x+phi) \rho(x-y) g(y) / Length
double Kernel::scalar(Fourier *f1, Fourier *f2, Complex *c, bool complem, int delta){ //phi -- shuffle phase
	if(f1->length !=length) return 0;
	if(f2->length !=length) return 0;
	double re=0, im=0;
	Fourier *zft=complem ? &cft : &ft;


	for(int i=0; i<length/2; i++){
		int idelta=(i+delta)%length;
		double RePhi=cos(2*PI*idelta/length);
		double ImPhi=sin(2*PI*idelta/length);
//		double RePhi=cos(2*PI*idelta*i/length);
//		double ImPhi=sin(2*PI*idelta*i/length);
		double ReF1 =f1->re[i] * RePhi - f1->im[i] *ImPhi;
		double ImF1 =f1->re[i] * ImPhi + f1->im[i] *RePhi;


		f1->re[i]=ReF1;
		f1->im[i]=ImF1;


		double RaRb_plus_IaIb =ReF1*f2->re[i] + ImF1*f2->im[i]; //==== Re(a)*Re(b)+Im(a)*Im(b)
		double RaIb_minus_IaRb=ReF1*f2->im[i] - ImF1*f2->re[i]; //==== Re(a)*Im(b)-Im(a)*Re(b)
		//==  Re=SUM Re(kern)(Re(a)*Re(b)+Im(a)*Im(b)) + Im(kern)*(Re(a)*Im(b)-Im(a)*Re(b))
		re+=(zft->re[i]*RaRb_plus_IaIb + zft->im[i]*RaIb_minus_IaRb)/length;
		//==  Im=SUM Im(kern)(Re(a)*Re(b)+Im(a)*Im(b)) - Re(kern)*(Re(a)*Im(b)-Im(a)*Re(b))
		im+=(zft->im[i]*RaRb_plus_IaIb - zft->re[i]*RaIb_minus_IaRb)/length;


//		int id=(i+delta)%length;
//		double ReF1 = f1->re[id];
//		double ImF1 = f1->im[id];
//
//		double RaRb_plus_IaIb =ReF1*f2->re[i] + ImF1*f2->im[i]; //==== Re(a)*Re(b)+Im(a)*Im(b)
//		double RaIb_minus_IaRb=ReF1*f2->im[i] - ImF1*f2->re[i]; //==== Re(a)*Im(b)-Im(a)*Re(b)
//		//==  Re=SUM Re(kern)(Re(a)*Re(b)+Im(a)*Im(b)) + Im(kern)*(Re(a)*Im(b)-Im(a)*Re(b))
//		re+=(zft->re[i]*RaRb_plus_IaIb + zft->im[i]*RaIb_minus_IaRb)/length;
//		//==  Im=SUM Im(kern)(Re(a)*Re(b)+Im(a)*Im(b)) - Re(kern)*(Re(a)*Im(b)-Im(a)*Re(b))
//		im+=(zft->im[i]*RaRb_plus_IaIb - zft->re[i]*RaIb_minus_IaRb)/length;
	}
	if(c!=0) {c->re=re; c->im=im;}


//	return sqrt(re*re+im*im);


	return re;
}


//======================================== Calculate distance (correlation)
double Kernel::dist(bool complem, int delta){
		return dist(&fx,&fy, complem, delta);
}




//======================================== Calculate distance (correlation)
double Kernel::dist(Fourier *f1, Fourier *f2, bool complem, int delta){
	Complex c0=Complex(), c1=Complex(), c2=Complex();


	double d12=scalar(f1,f2, &c0, complem, delta);
	double d11=scalar(f1,f1, &c1, true);
	double d22=scalar(f2,f2, &c2, true);
	Fourier *zft=complem ? &cft : &ft;


	dx11=d11; dx12=d12; dx22=d22; ex1=f1->re[0]; ex2=f2->re[0];
	//================================================ we should subtract the means
	double dd11=f1->re[0]*f1->re[0]*zft->re[0]/length;
	double dd12=f1->re[0]*f2->re[0]*zft->re[0]/length;
	double dd22=f2->re[0]*f2->re[0]*zft->re[0]/length;
	double e1=f1->re[0]/length, e2=f2->re[0]/length;
	//================================================ cumulative integral
	prod11+=d11;
	prod12+=d12;
	prod22+=d22;
	eprod1+=e1;
	eprod2+=e2;
	nprod++;
	//================================================ correlation for the window
	d11-=dd11; d12-=dd12; d22-=dd22;
	if(d11==0) return -400;
	if(d22==0) return -500;
	double cc=d12/sqrt(d11*d22);


	return cc;
}


//=========== inhibit Zero position of the kernel
double Kernel::NSCorrection(double x, double val){
	double x0=x*binSize;
	int zz=100/binSize; if(zz==0) zz=1;
	if(x0 < zz && x0>-zz){
		val*=(1-kernelNS);
	}
	return val;
}


//================================ General Kernel initiation
void Kernel::makeKernel(int n){
	init(n);
	double d=0;
	for(int i=0; i<length; i++){
		int x1=i, x2=i-length;
		ckern[length-i-1]=kern[i]=kernVal(x1)+kernVal(x2);
		d+=kern[i];
	}
	for(int i=0; i<length; i++){ckern[i]/=d; kern[i]/=d;}
	fft();
}


//======================== Normal Kernel ==============================
NormKernel::NormKernel():Kernel(){sigma=1; e=0; name= strdup("Normal_Kernel");}
NormKernel::NormKernel(double ee,double sgm, int l):Kernel(){
	sigma=sgm; e=ee; hasCompl=(e!=0);
	makeKernel(l);name= strdup("Normal_Kernel");
}


double NormKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma;
	double val=exp(-x*x/2);
	return NSCorrection(x0,val);
}


//======================== Left exp  Kernel ==============================
LeftExpKernel::LeftExpKernel():Kernel(){sigma=1; e=0; name= strdup("Left_Exp_Kernel");}
LeftExpKernel::LeftExpKernel(double ee,double sgm, int l):Kernel(){
	sigma=sgm; e=ee;
	makeKernel(l); name= strdup("Left_Exp_Kernel");
}
double LeftExpKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma; hasCompl=true;
	return NSCorrection(x0,(x>0)?0 : exp(x));
}




//======================== Right exp  Kernel ==============================
RightExpKernel::RightExpKernel():Kernel(){sigma=1;e=0;name= strdup("Right_Exp_Kernel");}
RightExpKernel::RightExpKernel(double ee,double sgm, int l):Kernel(){
	sigma=sgm; e=ee;
	makeKernel(l);name= strdup("Right_Exp_Kernel");
}
double RightExpKernel::kernVal(double x){
	double x0=x;
	x-=e; x/=sigma; hasCompl=true;
	return NSCorrection(x0,(x<0)?0:exp(-x));
}
//======================== Custom  Kernel ==============================
CustKernel::CustKernel():Kernel(){ initCust(0,1000);}
CustKernel::CustKernel(double ee,double sgm, int l):Kernel(){
	initCust(ee,sgm);
	makeKernel(l);name= strdup("custm_Kernel");
}
void CustKernel::initCust(double ee, double sgm){
	frml=frmlInit(customKern);
	frmlSetValue(frml,"sigma",sgm);
	frmlSetValue(frml,"e",ee);
	name= strdup("custm_Kernel");
}
double CustKernel::kernVal(double x){
	return NSCorrection(x,frmlCalc(frml,x));
}
//=============================================================================
//====================          Cross correlation              ================
//=============================================================================
//============================== initiation ===================================
void XYCorrelation::initXY(){
	init(kern->length);
	fx=&kern->fx; fy=&kern->fy;
	nCorr=nPlus=nMinus=0; min=max=av=sd=0;
	getMem0(correlation,length, "init correlation #1"); zeroMem(correlation,length);
	getMem0(corrMinus  ,length, "init correlation #2");	zeroMem(corrMinus  ,length);
	getMem0(corrPlus   ,length, "init correlation #3");	zeroMem(corrPlus   ,length);
	getMem0(datRe      ,length, "init correlation #4");	zeroMem(datRe      ,length);
	spectrumX=0; spectrumY=0;
}


//============================= Calculate XYCorrelation ======================
void XYCorrelation::calcXYCorr(int pos, bool cmpl1, bool cmpl2,  double cc){
	if(cc < -10) return;
	calcXYCorr(cmpl1, cmpl2);
	storeByChrom(pos, cc);
}


void XYCorrelation::calcXYCorr(bool cmpl1, bool cmpl2){
	for(int i=0; i<length; i++){
		double  ReX=fx->re[i],
				ReY=fy->re[i],
				ImX=fx->im[i],
				ImY=fy->im[i];
		ImX=cmpl1  ? -ImX : ImX;
		ImY=cmpl2  ? -ImY : ImY;
		datRe[i]=(ReX*ReY+ImX*ImY);
		datIm[i]=(-ReX*ImY-ImX*ReY);
//		datIm[i]=(-ReX*ImY+ImX*ReY);
	}
	calc(0);		//=========	reverse transformation
	norm();			//========= divide by length
	//========= normalize by std dev
}


//========== Store the cross-correlation by chromosomes
void XYCorrelation::storeByChrom(int pos, double corr){
	double delta=0.;
	Chromosome* chr=0;
	if(pos>=0) chr=getChromByPos(pos);
	double e1=track1->avWindow, d1=track1->sdWindow;
	double e2=track2->avWindow, d2=track2->sdWindow;


	for(int i=0; i<length; i++){
		double x=(re[i]-length*e1*e2)/(d1*d2*length);
		correlation[i]+=x;
		if(chr) {
			chr->distDens[i]+=x;
			if(corr >  delta) corrPlus [i]+=x;
			if(corr < -delta) corrMinus[i]+=x;
		}
	}
	if(chr){
		chr->densCount++;
		if(corr >  delta) nPlus++ ;
		if(corr < -delta) nMinus++;
	}
	nCorr++;
}


//================== Normalize the XY cross-correlation ====================
void XYCorrelation::normilize(){
	min=1.e+18; max=-1.e+18; av=sd=0;
	int n=0;
	for(int i=0; i<length; i++){
		double x=correlation[i]/nCorr;
		correlation[i]=x;
		if(nPlus) corrPlus [i]=corrPlus [i]/nPlus;
		if(nMinus) corrMinus[i]=corrMinus[i]/nMinus;


		if(min > x) min=x;
		if(max < x) max=x;
		av+=x; sd+=x*x; n++;
	}
	av/=n; sd=sd-av*av*n; sd/=n-1; sd=sqrt(sd);
}


void normChromDist(){
	for(int i=0; i<profWithFlanksLength; i++){
		for(int ich=0; ich<n_chrom; ich++){
			chrom_list[ich].distDens[i]/=chrom_list[ich].densCount;
		}
	}
}
void XYCorrelation::makeSpectrum()	{
	getMem(spectrumX, length,"spectX");
	getMem(spectrumY, length,"spectY");
	float *spX=fx->getSpectrum();
	float *spY=fy->getSpectrum();
	for(int i=0; i<length; i++){
		spectrumX[i] += spX[i];
		spectrumY[i] += spY[i];
	}
}













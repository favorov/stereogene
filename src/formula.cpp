/*
 * formula.cpp
 *
 *  Created on: 24 Jan 2017
 *      Author: Mironov
 */
#include "track_util.h"
#include <math.h>
struct Formula;
struct Identifier;
struct TrackNode;
//===================================================================
const int MULT	=0x1;
const int DIV	=0x2;
const int PLUS	=0x3;
const int MINUS	=0x4;
const int EQ	=0x5;
const int UMINUS=0xf;
const int SIN	=0x10;
const int COS	=0x11;
const int LOG	=0x12;
const int EXP	=0x13;
const int TAN	=0x14;
const int SQR	=0x15;
const int ABS	=0x16;
const int SIGN	=0x17;
const int ATAN	=0x18;


const int TRACK =0x30;	  	// node is track


const int CONST =0x80;	  	// node is a  constant
const int IDENT =0x81;	  	// node is an identifier


const int FunBeg  =0x81;  	// codes for functions
const int FunEnd  =0x100;  	// codes for functions
const int TrackBeg=0x101;  	// codes for tracks
const int TrackEnd=0x201;  	// codes for tracks
const int NodeBeg=0x1001;		// shift if nodes name
const int FunBBB =0x10;			// used in encoding b-formula


struct BFormula{
	static const char *functions[];
	static const char *shortFun;


	short bf[1024];
	int len;
	Formula *form;
	int root;


	BFormula(){init(0);}
	BFormula(Formula *ff){init(ff); }
	BFormula(BFormula *bf, int from, int to);


	void init(Formula *ff){len=0; form=ff;}
	void parseTerms(const char* s);
	const char *readIdent(const char *s, char*inent);
	const char *readNum(const char *s, char*num);
	const char *readTrack(const char *s, char*num);
	short getFun(const char *ident);
	bool isFun(int pos);
	bool isTrack(int pos);
	bool isNode(int pos);
	void replace(int from, int to, int val);
	int findBrace(int from);
	int parse();
	int parse(int from, int to);
	void parseOper(int pos, int oper);
	char * print(char*b);
	int getNode(int pos) 			{return bf[pos]-NodeBeg;}
	int getTrack(int pos) 		{return bf[pos]-TrackBeg;}
	int codeNode(int id)			{return id+NodeBeg;}
	void setNode(int pos, int id) 	{bf[pos]=codeNode(id);}
	int getFunct(int pos) 			{return bf[pos]-FunBeg+FunBBB;}
};
const char *BFormula::functions[]={"sin","cos","log","exp","tan","sqrt","abs","sign","atan",0};


const char *BFormula::shortFun="_sctelqag";
//===================================================================
//===================================================================


//===================================================================
//===================================================================
//===================================================================
//===================================================================
//===================================================================
Identifier::Identifier(Formula *frm, const char * idd){
	strcpy(name,idd);
	FNode *fn=new FNode(frm); fn->operation=IDENT;
	nodeID=fn->id;
	frm->identifiers[frm->nIds++]=this;
}
//===================================================================
TrackNode::TrackNode(Formula *frm, const char * trackName){
	strcpy(name,trackName);
	btr=new bTrack();
	btr->name=name;
	nodeID=frm->nTracks+TrackBeg;
	pos=-1; value=NA;
	frm->tracks[frm->nTracks++]=this;
}


double TrackNode::getValue(int p){
	if(p==pos) return value;
	pos=p;
	value=btr->getValue(pos,0);
	return value;
}


//===================================================================
char *TrackNode::print(char *b){
	sprintf(b,"[%s]",name);
	return b;
}
TrackNode::~TrackNode(){
//	if(btr) delete btr;
}
//===================================================================
void Formula::init(){
	formula=0;
	mainRoot=nNodes=nroots=nTracks=nIds=0;
	arg   = fNodes[addIdent("x")];
	e  	  = fNodes[addIdent("e")];
	sigma = fNodes[addIdent("sigma")];
}


Formula::~Formula(){
	if(formula!=0) {free(formula);} formula=0;
	for(int i=0; i<nTracks; i++){
		if(tracks[i]) del(tracks[i]);
		tracks[i]=0;
	}
}
//===================================================================
double Formula::getValue(const char* ident){
	Identifier *idx=getIdentificator(ident);
	if(idx==0) return NAN;
	return fNodes[idx->nodeID]->value;
}
//===================================================================
void Formula::setArg(double x){
	arg->value=x;
	for(int i=0; i<nTracks; i++){
		int pp=(int)x;
		tracks[i]->getValue(pp);
	}
}


//===================================================================
void Formula::setValue(const char* ident, double v){
	Identifier *idx=getIdentificator(ident);
	int fnId= (idx==0) ? addIdent(ident) : idx->nodeID;
	fNodes[fnId]->value=v;
}


//===================================================================
Identifier *Formula::getIdentificator(const char *ident){
	for(int i=0; i< nIds; i++)
		if(strcmp(identifiers[i]->name,ident)==0) return identifiers[i];
	return 0;
}
//===================================================================
TrackNode *Formula::getTrack(const char *ident){
	for(int i=0; i< nTracks; i++)
		if(strcmp(tracks[i]->name,ident)==0) return tracks[i];
	return 0;
}
//===================================================================
int Formula::addIdent(const char* ident){
	Identifier *idx=getIdentificator(ident);
	if(idx==0)  return (new Identifier(this, ident))->nodeID;
	else 	    return idx->nodeID;
}
//===================================================================
int Formula::addTrack(char* ident){
	ident=trim(ident);
	TrackNode *idx=getTrack(ident);
	if(idx==0)  return (new TrackNode(this, ident))->nodeID;
	else 	    return idx->nodeID;
}
//===================================================================
void Formula::parse(const char* input){
	char b[1024], *s;
	strcpy(b,input);
	formula=strdup(input);
	s=strtok(b,";");
	while(s){
		BFormula bf(this);
		bf.parseTerms(s);
		int r=bf.parse();
		roots[nroots++]=r;
		if(strchr(s,'=')==0) mainRoot=r;
		s=strtok(0,";");
	}
}
//===================================================================
double Formula::calc(double x){
	setArg(x);
	double y=calc();
	return y;
}
//===================================================================
double Formula::calc(){


	for(int i=0; i<nroots; i++) {
		getNode(roots[i])->calc();
	}
	return  getNode(mainRoot)->value;
}
//===================================================================
//===================================================================
//===================================================================
FNode::FNode(Formula* form){
	formula=form;
	childL=childR=-1;
	operation=0;
	value=NAN;
	id=form->nNodes++;
	form->fNodes[id]=this;
}
//===================================================================
char *FNode::print(char *b){
	const char *oprs="*/+-=_sctelqa:u@#~!";
	const int opcodes[]=
		{MULT,DIV,PLUS,MINUS,EQ,UMINUS,SIN,COS,TAN,EXP,LOG,SQR,ABS,SIGN,ATAN,IDENT,CONST,TRACK};
	char oper='?';
	for(int i=0; oprs[i]; i++) if(operation==opcodes[i]) oper=oprs[i];
	char val[20], chL[10], chR[10];
	if(childL <0) strcpy(chL," "); else snprintf(chL,sizeof(chL),"%i",childL);
	if(childR <0) strcpy(chR," "); else snprintf(chR,sizeof(chR),"%i",childR);
	sprintf(val,"%f",value);
	sprintf(b,"$%i\tv=%s\t$%s\t%c\t$%s %x",id, val, chL,oper,chR,operation);
	if(operation==IDENT){
		for(int i=0; i<formula->nIds; i++){
			if(formula->identifiers[i]->nodeID==id)
				{sprintf(b+strlen(b),"\t%s",formula->identifiers[i]->name); break;}
		}
	}
	if(operation==TRACK){
		for(int i=0; i<formula->nTracks; i++){
			if(formula->tracks[i]->nodeID==id)
				sprintf(b+strlen(b),"\t[%s]",formula->tracks[i]->name);
		}
	}
	return b;
}
//===================================================================
double FNode::calc(){
	double v1=NAN, v2=NAN;
	FNode *fn1=0, *fn2=0;
	if(childL >=0){fn1=formula->getNode(childL); v1=fn1->calc();}
	if(operation != TRACK)
	if(childR >=0){fn2=formula->getNode(childR); v2=fn2->calc();}
	switch(operation){
	case MULT 	:  value=v1*v2; 	break;
	case DIV 	:  value=v1/v2; 	break;
	case PLUS 	:  value=v1+v2; 	break;
	case MINUS 	:  value=v1-v2; 	break;
	case EQ 	:  value=fn1->value=fn2->value; break;
	case UMINUS :  value=-v1;		break;
	case SIN 	:  value=sin(v1);	break;
	case COS 	:  value=cos(v1);	break;
	case LOG 	:  value=log(v1);	break;
	case EXP 	:  value=exp(v1);	break;
	case TAN 	:  value=tan(v1);	break;
	case SQR 	:  value=sqrt(v1);	break;
	case ABS 	:  value=abs(v1);	break;
	case SIGN 	:  value=sign(v1);	break;
	case ATAN 	:  value=atan(v1);	break;
	case TRACK	:
			value=formula->tracks[childR]->getValue((int)v1); 	// in the track childR=track no
		 break;
	}
	return value;
}
//===================================================================
//===================================================================
//===================================================================
const char *BFormula::readIdent(const char *s, char*ident){//=== read identifier
	for(char *ss=ident; *s && (*s=='_' || isalnum(*s));) {*ss++=*s++; *ss=0;}
	return s;
}


const char *BFormula::readTrack(const char *s, char*ident){//=== read track
	s++;
	for(char *ss=ident; *s && *s!=']';) {*ss++=*s++; *ss=0;}
	if(*s) s++;
	return s;
}
//===================================================================
const char *BFormula::readNum(const char *s, char*num){//==== read a number
	char *t=num;
	for(;*s!=0 && isdigit(*s);) {*t++=*s++;} *t=0;
	if(*s=='.')              	{*t++=*s++;} *t=0;
	for(;*s!=0 && isdigit(*s);) {*t++=*s++;} *t=0;
	if(*s=='e' || *s=='E'){  	{*t++=*s++;} *t=0;
		if(*s=='+' || *s=='-')	{*t++=*s++;} *t=0;
		for(;*s!=0 && isdigit(*s);){*t++=*s++;} *t=0;
	}
	return s;
}
//===================================================================
short BFormula::getFun(const char *ident){//======= get a function from the list
	for(int i=0; functions[i]; i++){
		if(strcmp(ident,functions[i])==0) return i+FunBeg;
	}
	return 0;
}
//===================================================================
void BFormula::parseTerms(const char* input){//===== convert text formula to a binary formula
	char ident[256];
	len=0;
	for(const char *s=input; *s;){
		if(*s=='_' || isalpha(*s)){	//================ start identifier
			s=readIdent(s,ident);
			short fun=getFun(ident);
			if(fun) bf[len++]=fun;
			else	setNode(len++,form->addIdent(ident));
		}
		else if(*s=='['){			//=============== track
			s=readTrack(s,ident);
			bf[len++]=form->addTrack(ident);
			if(*s!='('){
				bf[len++]='(';
				setNode(len++,form->arg->id);
				bf[len++]=')';
			}
		}
		else if(*s=='.' || isdigit(*s)){
			s=readNum(s,ident);
			FNode *fn=new FNode(form); fn->value=atof(ident); fn->operation=CONST;
			setNode(len++,fn->id);
		}
		else if(isspace(*s)) s++;
		else bf[len++]=*s++;
	}
}
//===================================================================
char * BFormula::print(char *bb){//======= print binary formula to a buffer
	char *b=bb;	*b=0;
	char bx[40];
	sprintf(bb,"len=%i  ",len); b=bb+strlen(bb);
	for(int i=0; i<len; i++){
		if(isNode(i)) 		{sprintf(bx,"$%x",bf[i]); strcat(b,bx); b+=strlen(b);}
		else if(isTrack(i)) {sprintf(bx,"^%x",bf[i]); strcat(b,bx); b+=strlen(b);}
		else if(isFun(i)) 	{*b++=shortFun[bf[i]-FunBeg+1]; *b=0;}
		else      			{*b++=bf[i]; *b=0;}
	}
	return bb;
}
//====================================================================
void BFormula::replace(int from, int to, int val){//== replace fragment by a node
	for(int i=to, j=from+1; i<len && j<len; i++,j++) bf[j]=bf[i];
	int ll=to-from-1; len-=ll;
	bf[from]=val;
}
BFormula::BFormula(BFormula *bfx, int from, int to){//======== make fragment of given binary formula
	init(bfx->form);
	len=to-from;
	for(int i=from, j=0; i < to; i++, j++) bf[j]=bfx->bf[i];
}
//=====================================================================
int BFormula::findBrace(int from){//========== find appropriate brace
	int level=0;
	for(int i=from; i<len; i++){
		if(bf[i]=='(') level++;
		if(bf[i]==')') {level--; if(level==0) return i;}
	}
	errorExit("Formula syntax error #3 (brackets eror): \"%s\"",form->formula);
	return -1;
}
//=====================================================================
bool BFormula::isFun(int pos){//================= is the element function
	return (bf[pos] >= FunBeg) && (bf[pos] <= FunEnd);
}
//=====================================================================
bool BFormula::isTrack(int pos){//================= is the element function
	return (bf[pos] >= TrackBeg) && (bf[pos] < TrackEnd);
}


//=====================================================================
bool BFormula::isNode(int pos){//================= is the element Node
	return bf[pos] >= NodeBeg;
}
//=====================================================================
void BFormula::parseOper(int pos, int oper){//=========== parse operation
	FNode *fn=new FNode(form);
	fn->childL=getNode(pos-1);
	fn->childR=getNode(pos+1);
	fn->operation=oper;
	if(oper==EQ){
		FNode *left=form->getNode(fn->childL);
		if(left->operation!=IDENT) 	errorExit("Formula syntax error #4 (left to \'=\' is not identifier): \"%s\"",form->formula);


	}
	replace(pos-1,pos+2,codeNode(fn->id));
}
//=====================================================================
int BFormula::parse(int from, int to){//============ parse fragment. return: the node
	BFormula bf1(this,from, to);
	bf1.parse();
	return bf1.bf[0];
}
//=====================================================================
int BFormula::parse(){
	int from=0, to=len;
	while(bf[from]=='(' && bf[to-1]==')') {from++; to--;}
	for(int i=0, j=from; j<to; i++,j++) {bf[i]=bf[j];} len=to-from;
	for(int i=0; i < len; i++){//========== find brakets
		if(bf[i]=='('){
			int j=findBrace(i);
			int node= (j-i==2) ? bf[i+1] : parse(i+1,j);
			replace(i,j+1, node);
		}
	}
	for(int i=0; i < len-1; i++){//======== parse functions
		if(isFun(i)){
			if(!isNode(i+1)) errorExit("synax error #2 (unknown function): \"%s\"",form->formula);
			FNode *fn=new FNode(form); fn->childL=getNode(i+1);
			fn->operation=getFunct(i);
			replace(i,i+2,codeNode(fn->id));
		}
	}


	for(int i=0; i < len-1; i++){//======== parse tracks
		if(isTrack(i)){
			if(!isNode(i+1)) errorExit("synax error #2 (track without arg): \"%s\"",form->formula);
			FNode *fn=new FNode(form); fn->childL=getNode(i+1);
			fn->childR=getTrack(i);
			fn->operation=TRACK;
			replace(i,i+2,codeNode(fn->id));
		}
	}
	for(int i=0; i < len-1; i++){//========= parse unitar '-'
		if(bf[i]=='-' && (i==0 || !isNode(i-1))){
		FNode *fn=new FNode(form); fn->childL=getNode(i+1);
		fn->operation=UMINUS;
		replace(i,i+2,codeNode(fn->id));
		}
	}
	for(int i=1; i < len-1; i++){//=========== parse *, /
		if(bf[i] == '*') {parseOper(i,MULT ); i--;}
		if(bf[i] == '/') {parseOper(i,DIV  ); i--;}
	}
	for(int i=1; i < len-1; i++){//=========== parse +-
		if(bf[i] == '+') {parseOper(i,PLUS ); i--;}
		if(bf[i] == '-') {parseOper(i,MINUS); i--;}
	}
	for(int i=1; i < len-1; i++){//=========== parse =
		if(bf[i] == '=') parseOper(i,EQ   );
	}
	if(len > 1) {errorExit("syntax error #1: \"%s\"",form->formula);}
	return getNode(0);
}
//=====================================================================
//=====================================================================
//===================================================================== frmlInit
//=====================================================================
Formula *frmlInit(const char* txt){
	Formula *f=new Formula();
	f->parse(txt);
	return f;
}
void 	frmlClose	(Formula* f){del(f);}
double 	frmlCalc	(Formula* f, double x){return f->calc(x);}
void 	frmlSetValue(Formula* f, const char* txt, double val){f->setValue(txt,val);}
double 	frmlGetValue(Formula* f, const char* txt){return f->getValue(txt);}




void test_formula(){
	const char *input="a=2; [h3k4](x+3)+a";
//	const char *input="s1=0.5; s2=1; m1=1; m2=-1; y1=(x-m1)/s1; y2=(x-m2)/s2; 1/s1*exp(-y1*y1) + 1/s2*exp(-y2*y2)";


	Formula formula;


	formula.parse(input);


//char b[256];
//for(int i=0; i<formula.nNodes; i++){
//	deb(formula.getNode(i)->print(b));
//}


double x=80;
double y=formula.calc(x);
printf("f(%f)=%f\n",x,y);
//deb("f(%f)=%f",x,y);


//for(int i=0; i<100; i++){
//	deb("%i\t%f",i,formula.calc(i));
//}




//	formula.setValue("sigma",0.2);
//	formula.setValue("x",0);
//	formula.setValue("y",0);
//
//	double d=0.00001, x=0, y=0;
//	Timer t;
//	for(int i=0; i<10000000; i++){
//		y=y+exp(-x)*d; x=x+d;
//	}
//	deb("%s",t.getTime());
//	t.reset();
//	for(int i=0; i<10000000; i++){
//		formula.calc();
//	}
//	deb("%s",t.getTime());
//	deb("%f\t%f",x,y);
//	deb("%f\t%f",formula.getValue("x"),formula.getValue("y"));
//	for(double x=-3; x<3; x+=0.1){
//		deb("%f\t%f",x, formula.calc(x));
//	}
exit(0);
}



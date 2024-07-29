#include "track_util.h"


//=========================================
//=========================================
//=========================================
//==========================================================================================



void testMtx(){
	int n=3;
	double mm[]={	1,   -0.2,	0.3,
				 -0.2,    1,	0.1,
				  0.3,	0.1,	1
	};
	VectorX *v=new VectorX(n);
	Matrix *m=new Matrix(n,mm);
	m->printMtx();

	m->eigen(v);
	v->print();

	exit(0);
}


const char * progName="Confounder";
const int progType=SG;


void printProgDescr(){
	printf("\n");
	printf("The Confounder program creates a confounder track using set of the tracks\n");
	printf("Usage:\n");
	printf("$ ./parse_genes [-parameters] [RefSeq or GENECODE file]\n");
	printf("\n");
}


void printMiniHelp(){
	printf("\n");
	printf("The Confounder program creates a confounder track using set of the tracks\n");
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./%s [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n",progName);
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}


int main(int argc, char **argv) {
	debugFg=3; if(debugFg) clearDeb();
//	testMtx();


	initSG(argc, argv);


	Preparator();
	Covariator();
	fflush(stdout);
	fclose(stdout);
	return 0;
}

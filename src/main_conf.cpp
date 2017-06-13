#include "track_util.h"

const char * progName="StereoGene";
const int progType=SG;


void printMiniHelp(){
	printf("\n");
	printf("The %s program compares pairs of tracks and calculates kernel correlations\n",progName);
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./%s [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n",progName);
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}

int main(int argc, char **argv) {
	initSG(argc, argv);
	if(debugFg) clearDeb();

	Preparator();
	if(confFile==0) confFile=strdup("confounder");
	Covariator();
	fflush(stdout);
	fclose(stdout);
	return 0;
}

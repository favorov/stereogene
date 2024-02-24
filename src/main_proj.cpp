/*
 * main_proj.cpp
 *
 *  Created on: 03 Jan 2017
 *      Author: andrey
 */
#include "track_util.h"


const char * progName="Projection";
const int progType=SG;




void printProgDescr(){
	printf("\n");
	printf("The Projection program creates tracks with exclusion of the defined confounder\n");
	printf("Usage:\n");
	printf("$ ./Projection [-parameters] track1 track2 ...\n");
	printf("\n");
}
void printMiniHelp(){
	printf("\n");
	printf("The Projection program creates tracks with exclusion of the defined confounder\n");
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./Projection [-parameters] track1 track2 ...\n");
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}


int main(int argc, char **argv) {
	initSG(argc, argv);
	if(debugFg) {clearDeb(); debugFg=DEBUG_LOG|DEBUG_PRINT;}


	Preparator();


	Projector();
	fflush(stdout);
	fclose(stdout);
	return 0;
}









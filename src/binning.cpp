/*
 * main_proj.cpp
 *
 *  Created on: 03 Jan. 2017
 *      Author: Mironov
 */
#include "track_util.h"


const char * progName="Binning";
const int progType=SG;




void printProgDescr(){
	printf("\n");
	printf("The Binning program create a track with binned data\n");
	printf("Usage:\n");
	printf("$ ./binning [-parameters] track1 track2 ...\n");
	printf("\n");
}
void printMiniHelp(){
	printf("\n");
	printf("The Binning program create a track with binned data\n");
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./binning [-parameters] track1 track2 ...\n");
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}








void binning(const char *fname){
	bTrack *tmp=new bTrack();
	tmp->writeBinnedProf(fname);


	del(tmp);
	if(fProfile) del(fProfile);
}


int main(int argc, char **argv) {
	initSG(argc, argv);
	if(debugFg) {clearDeb(); debugFg=DEBUG_LOG|DEBUG_PRINT;}


	for(int i=0; i<nfiles; i++){
		char *fname=files[i].fname;
		if(fname==0 || strlen(trim(fname))==0) continue;
		binning(fname);
	}


	fflush(stdout);
	fclose(stdout);
	return 0;
}









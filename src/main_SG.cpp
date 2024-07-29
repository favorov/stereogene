#include "track_util.h"

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

//==============================================================================
// for debugging :
// set debugFg=DEBUG_LOG|DEBUG_PRINT
// set debS string for module identification//
// Use:  deb(n);  // print debug information as number
//       deb(format,...)    // print debug information as printf
//       deb(n,format,...)  // print debug information as a number and printf
// example:
// debS="fun1";
// deb(1);
// ....
// deb(2,"%i %f", n, d);
// ....
// deb("OK");


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


void printProgDescr(){
	printf("\n");
	printf("The StereoGene program compares pairs of tracks and calculates kernel correlations\n");
	printf("Usage:\n");
	printf("$ ./StereoGene [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n");
	printf("\n");
}


//============================================ Tests =========================================
void test(const char* path){
	clearDeb();
	debugFg=DEBUG_LOG|DEBUG_PRINT;

	readInt("99");
	readInt("99kb");
	readInt("99 kb");
	readInt("99 kb");
	readInt("99 k");
	readInt("99Mb");

exit(0);
}

//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================

int main(int argc, char **argv) {
//	test("../Tracks/RCDB_all_to_all");
	clearDeb();	debugFg=DEBUG_LOG|DEBUG_PRINT;

	initSG(argc, argv);
	writeLog("====== Start ====== deb=%i\n",debugFg);

//===========================================
	Preparator();
	Correlator();
	fflush(stdout);
	writeLog("====== End ======\n");
	return 0;
}

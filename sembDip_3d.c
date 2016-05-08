#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defns.h"
#include "segy.h"
#include "segyIO_class.h"
#include <time.h>
#include <omp.h>

#define FLUSH while (getchar() !='\n');

/***********************************************************************
 * Function: main
 * Input/output:  input seismic file name, output dip file name, 
 * 		output dip file name, number of traces used in semblance calc,
 * 		dip as max ms/trace allowable, semblance window size
 * 	      
 * Descrip: Compute dips via semblance of adjacent traces.  
 * 		This dip is output into separate files to be read in at
 * 		a later time, and (hopefully) used to smear the velocities
 * 		along dips/structure. 
 * ********************************************************************/
 
void sembl(float ***trGth, int ns, int numTr, int iNdx, int window, int maxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx,  float thresh, int xMax, int rate, int reso);

int main(int argc, char *argv[])
{
FILE  *fin=NULL,  *foutX=NULL, *foutI=NULL, *fsemb=NULL; //fout = file out , fmodel = model input
int   fmt, i, itr, ierr, k;
int   ns;	/* Number of samples in a trace	*/
float dts;	/* Sample rate in seconds	*/
segy  thdr;	/* Segy trace header		*/
segy thdrPad;
bhed  bh;	/* Segy file binary header	*/
char chdr[3202];	/* segy file character header */
int   endian=1; /* Which side of the egg	*/
int   status=0; /* Return status for segy calls, I should check*/
char  txthdr[3201];	/* The character file header*/

int choice;
char buff[80];
int first =1; /*value to track whether this is first gther read in*/
float ***trGth; /*2D array for gather used in semblance-dip calc*/
segy ***trGthHdr;
float *tempTr; /*holding place for traces as they are read in*/
float ***dipTrI;
float ***dipTrX;
float ***semTr;
int *xNdxMx; /*Keeps track of number of xline on each inline*/
int *xNdxMn; /*first xline on each inline*/
int statusR=0; /*Status of last segyRead*/
int statusW=0; /*Status of last segyWrite*/
int numTr=3; /*Radius of traces, about the trace being examined, used in
             the semblance-dip calculation*/
int maxDip=6; /*maximum allowed dip, IN SAMPLES PER TRACE*/
int window=40; /*Size of sembalnce window in number of samps*/
int shift = 0; /*shifts starting sample point down for dip calc.  if you
				started at the 0th sample, the dip calc would run into
				* negative sample points*/
int iOrigin = 0; /*origin of inline grid*/
int xOrigin = 0; /*origin of xline grid*/
int iMax = 0;
int xMax = 0;
int line=0;
int iline = 0;
int xline = 0;
int j;
int ilineOld;			
int ilineCnt=0;
float thresh = 0;
int skip = 0;
int rate;
int reso;

int num10Percent;

float temp;
time_t now;
	time(&now);

endian = checkEndian();

thdr.iline = thdr.ep;
thdr.xline=thdr.cdpt;


	if(argc ==14){
/*Open fileS as readable/writable binary*/
		fin = fopen(argv[1], "rb");
		if (fin == NULL) {
			printf("Unable to open the Input file.  Please check the name.\n");
		return -1;
		}
		foutI=fopen(argv[2],"r");
		if(foutI!=NULL){
			printf("\ndipI output file Seems to exist. ");
			printf("\nWould you like to overwrite it (-1 = yes)?");
			fgets(buff,80,stdin);
			sscanf(buff,"%d",&choice);
			if (choice == -1){
				fclose(foutI);
				foutI=fopen(argv[2],"wb");
			}else{
				fclose(foutI);
				return 0;
			}
		//FLUSH;
		choice=0;
		}else {
			foutI=fopen(argv[2],"wb");
		}
		
		foutX=fopen(argv[3],"r");
		if(foutX!=NULL){
			printf("\ndipX output file Seems to exist. ");
			printf("\nWould you like to overwrite it (-1 = yes)?");
			fgets(buff,80,stdin);
			sscanf(buff,"%d",&choice);
			if (choice == -1){
				fclose(foutX);
				foutX=fopen(argv[3],"wb");
			}else{
				fclose(foutX);
				return 0;
			}
		//FLUSH;
		choice=0;
		}else {
			foutX=fopen(argv[3],"wb");
		}	
			
		fsemb=fopen(argv[4],"r");
		if(fsemb!=NULL){
			printf("\nsemblance output file Seems to exist. ");
			printf("\nWould you like to overwrite it (-1 = yes)?");
			fgets(buff,80,stdin);
			sscanf(buff,"%d",&choice);
			if (choice == -1){
				fclose(fsemb);
				fsemb=fopen(argv[4],"wb");
			}else{
				fclose(fsemb);
				return 0;
			}
		//FLUSH;
		choice=0;
		}else {
			fsemb=fopen(argv[4],"wb");
		}

		if(sscanf(argv[5],"%d", &numTr)!=1){
			printf("\nnumTr entered not valid. Exiting... ");
			return -1;
		}
		if(sscanf(argv[6],"%d", &maxDip)!=1){
			printf("\nmaxDip entered not valid. Exiting... ");
			return -1;
		}
		if(sscanf(argv[7],"%d", &window) !=1){
			printf("\nwindow entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[8],"%f", &thresh) !=1){
			printf("\nwindow entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[9],"%d", &reso) !=1){
			printf("\nwindow entered not valid. Exiting...");
			return -1;
		}
		
		if(sscanf(argv[10], "%d", &iOrigin) != 1){
			printf("\niOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[11], "%d", &xOrigin) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[12], "%d", &iMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[13], "%d", &xMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}		
	} else {
	fprintf(stderr,"\n***************************************************************************\n");
      fprintf(stderr,"Computes the dip and coherancy at each sample via semblance. \n\n");
      fprintf(stderr,"Program expects the following command line: \n ");
      fprintf(stderr,"sembDip <in.sgy><dipI.sgy><dipX.sgy><sem.sgy><Num><maxDip><win><thresh>\n");
      fprintf(stderr,"<reso><iOrigin><xOrigin><iMax><xMax>\n\n");
      fprintf(stderr,"in.sgy: Input filename. \n\n");
      fprintf(stderr,"dipI.sgy: Output Inline dip sgy filename.\n\n");
      fprintf(stderr,"dipX.sgy: Output Xline dip sgy filename.\n\n");
      fprintf(stderr,"sem.sgy: Output semblance  sgy filename.\n\n");
      fprintf(stderr,"Num: Number of traces for semblance given as  \n");
      fprintf(stderr,"     distance from trace being analysed. \n");
      fprintf(stderr,"     ex. numTr = 1 uses 9 traces. 3 by 3 block.\n");
      fprintf(stderr,"                 o--o--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--x--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--o--o\n\n");
      fprintf(stderr,"maxDip: Maximum dip in ms/trace.\n\n");
      fprintf(stderr,"win: Number of samps used for semblance analysis.\n\n");
      fprintf(stderr,"thresh: minimum value for semblance - dips zerored if below this value.\n\n");
      fprintf(stderr,"reso: resolution in points/samp. Ex. reso = 1 calculates dip to within one \n");
      fprintf(stderr,"      point per sample per trace.  reso = 2 calcs to within 2 points/samp.\n\n");
      fprintf(stderr,"iOrigin: Starting inline number.\n\n");
      fprintf(stderr,"xOrigin: Starting xnline number.\n\n");
      fprintf(stderr,"iMax: Max iline. \n\n");
      fprintf(stderr,"xMax: Max xline.");
    fprintf(stderr, "\n**************************************************************************\n");
      return 0;
    } 
    segyReadHeader(fin, chdr, &bh, endian);
	ns=bh.hns;
	rate = bh.hdt;
	rate=rate/1000;

	/*print parameters. ask user to verify they are correct-option to abort*/
	printf("\nInput paramters are:");
	printf("\n\n  Input seismic: %s", argv[1]);
	printf("\n  Output dipI: %s", argv[2]);
	printf("\n  Output dipX: %s", argv[3]);
	printf("\n  Output semb: %s", argv[4]);
	//printf("\n Output velocity:  %s", fvelout);
	printf("\n  Trace radius: %d", numTr);
	printf("\n  Max dip (in samples!): %d", maxDip);
	printf("\n  Window (in samples!): %d", window);
	printf("\n  threshhold: %f", thresh);
	printf("\n  resolution: %d", reso);
	printf("\n  inline origin: %d", iOrigin);
	printf("\n  xline origin: %d", xOrigin);
	printf("\n  iMax : %d", iMax);
	printf("\n  xMax: %d", xMax);
	printf("\n\n  With a samp rate of %d ms(ft)/trace and resolution of %d pts/samp, the dip", rate,reso);
	printf("\n  will be computed to within %2.1f ms(ft)/trace with a max dip of %d ms(ft).", (float)rate/reso,maxDip*rate);
	
	printf("\n\nAre these correct? (-1: exit)");
	fgets(buff,80,stdin);
	sscanf(buff,"%d",&choice);
	if(choice == -1){ return 0;} 

	if(segyWriteHeader(foutI, &chdr, &bh, endian)!=0){
		printf("\n Unable to write header for output data!! \n");
		return -1;
	}
	if(segyWriteHeader(foutX, &chdr, &bh, endian)!=0){
		printf("\n Unable to write header for output data!! \n");
		return -1;
	}
	if(segyWriteHeader(fsemb, &chdr, &bh, endian)!=0){
		printf("\n Unable to write header for output data!! \n");
		return -1;
	}	
	
	/*Converting numTr to the width of cube ex. 3x3*/
	numTr=(numTr)*2+1;
	/*I prefer the window to be an odd number so as to be symmetric*/
	if(window % 2 ==0){
		window = window + 1; 
	}
	/*allocating memory for the traces used for semblance calc/analysis*/
	trGth = (float***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         trGth[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 trGth[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	dipTrI = (float ***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         dipTrI[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 dipTrI[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 dipTrX = (float ***)calloc((numTr) , sizeof(float **));
     for (i=0; i<(numTr); ++i){
         dipTrX[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 dipTrX[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }	
	 semTr = (float ***)calloc((numTr) , sizeof(float **));
     for (i=0; i<(numTr); ++i){
         semTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 semTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }	
	 trGthHdr = (int ***)calloc((numTr) , sizeof(int **));
	 for (i=0; i<(numTr); ++i){
         trGthHdr[i] = (int **)calloc((xMax-xOrigin+2),sizeof(int*)); 
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 trGthHdr[i][j] = (int *) calloc((1),sizeof(segy));
		 }
	 }

	 tempTr = (float *)calloc((ns+1),sizeof(float));
	 xNdxMx = (int *)calloc ((iMax - iOrigin +2), sizeof(float)); //to track the starting and ending xlines
	 xNdxMn = (int *)calloc ((iMax - iOrigin +2), sizeof(float)); //to track the starting and ending xlines

	 /*Reading in frist Trace*/
	statusR = segyReadTrace(fin, &bh, &thdr,tempTr, ns, endian);
	
	/*im setting the iOrigin to the first inline read in-will be useful later*/
	if(thdr.iline != iOrigin){
		printf("\n\nFirst inline read does not match iOrigin: %d\n",iOrigin);
		printf("Setting iOrigin equal to first iline read in!: %d\n\n",thdr.iline);
		iOrigin=thdr.iline;
	}
	/*are the values sensible?*/
	if (thdr.iline > iMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}else if (thdr.xline - xOrigin < 0){
		printf("\n\nFirst xline read less than xOrigin.\n\n");
		return 0;
	}else if (thdr.xline > xMax){
		printf("\n\nFirst xline read great than xMax.\n\n");
		return 0;
	}
	
	/*save the xline values and set iline switch - these are the index numbers, so
	 * for example if the current thdr.xline =101 and xOrigin =101 then xNdxMn =0, for
	 * that particular inline*/
	xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	ilineOld=thdr.iline;
	
	/*read in the first numTr number of traces*/
	while ( thdr.iline - iOrigin < numTr){

		/*copy the values read from previous trace read*/
		memcpy(trGth[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTr,(ns+1)*sizeof(float));
		memcpy(trGthHdr[thdr.iline - iOrigin][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));

		statusR = segyReadTrace(fin, &bh, &thdr,tempTr, ns, endian);	
		
		/*check to make sure data read is good, otherwise abort*/
		if (thdr.iline > iMax){
			printf("\n\nFound iline greater than iMax.\n\n");
			return 0;
		}else if (thdr.iline < iOrigin){
			printf("\n\nFound iline less than iOrigin.\n\n");
			return 0;
		}else if (thdr.xline - xOrigin < 0){
			printf("\n\nFound xline  less than xOrigin.\n\n");
			return 0;
		}else if (thdr.xline > xMax){
			printf("\n\nFound xline greater than xMax.\n\n");
			return 0;
		}
			
		/*if it gets to end of reading before having read in enough traces abort
		*else end o file and you have enough lines exit while loop - overKill! */
		if(statusR!=0){
			if(thdr.iline - iOrigin < numTr -1 ){
				printf("\nNot enough ilines for a single semblance calc!\n");
				return  0;
			}else{	
				break; /*minimum number of lines met...exit loop*/
			}
		}
		
		    /*if the xline number just read is large than xNdxMx set xNdxMx to thdr.xline just read*/
			if ( (thdr.xline-xOrigin) > xNdxMx[thdr.iline-iOrigin]) xNdxMx[thdr.iline - iOrigin]=thdr.xline - xOrigin;
			/*make note of iline changes*/
			if (ilineOld != thdr.iline){
				xNdxMn[thdr.iline-iOrigin]=thdr.xline-xOrigin;
				ilineOld=thdr.iline;
			}
	}
	
	/*ilineCnt starts at the first iline being Processed (numTr/2 + iOrigin) and increments after
	 * each read. This could be made less confusing. */
	ilineCnt=numTr/2;
	
	printf("\n\n ****************");
			time(&now);
			printf("\n Start Time: %s Percent Done: %d",  ctime(&now),0 );

	/*compute semblance and dips could send xNdxMn[ilineCnt] rather than sending both ilineCnt and *xNdxMn */
	sembl(trGth, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,thresh,xMax,rate,reso);

	num10Percent = (iMax-iOrigin + 1)/10;
	if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
		printf("\n\n ****************");
		time(&now);
		printf("\n Finished Inline: %d Time: %s Percent Done: %d ", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline - iOrigin + 1)/num10Percent*10 );
	}

	/*spits out first numTr/2 +1  lines*/
	for (i=0;i<=numTr/2;++i){
		for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		//for (j = 0; j <= xMax-xOrigin; ++j){ //Use this for padded Output
		
			/*if no header exists (example padded trace within inline with iline = xline =0) put inline and xline numbers into headers*/
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
				segyWriteTrace(foutI, trGthHdr[i][j], dipTrI[numTr/2][j], ns, endian);
				segyWriteTrace(foutX, trGthHdr[i][j], dipTrX[numTr/2][j], ns, endian);
				segyWriteTrace(fsemb, trGthHdr[i][j], semTr[numTr/2][j], ns, endian);
			}
	}
	
	/*while the segyReadTrace function is returning good reads (statusR==0), continue reading*/
   while (statusR == 0){ 

		/*Shift all the traces back by one to make room for next trace, meaning that
		* trGth[1] -> trGth[0], trGth[2] -> [1], and tempTr -> trGth[numTr-1]*/
		for (k =0 ;k<numTr-1;++k){
			for(j=0;j<xMax-xOrigin+1;++j){
					memcpy(trGth[k][j],trGth[k+1][j],(ns+1)*sizeof(float));
					memcpy(trGthHdr[k][j],trGthHdr[k+1][j], sizeof(*trGthHdr[0][0]));
			}
		}
		
		/*zeroing out last inline of gather before reading new data into it - just in case!*/
		for (j =0; j < xMax-xOrigin +1;++j){
			memset(trGth[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(trGthHdr[numTr-1][j],0,sizeof(*trGthHdr[0][0]));
		}

		/*Copying the last read tempTr and and thdr before continuing with rest of inline. This is done since the previous
		 * inline read stops when it finds the first trace of the next inline, meaning that the tempTr and thdr still contain
		 * the first trace of this inline when the inline read function is called next*/
		while(statusR ==0 && thdr.iline ==ilineOld){
			
			/*copy the previously read trace values*/
			memcpy(trGth[numTr - 1][thdr.xline - xOrigin], tempTr, (ns + 1)*sizeof(float));
			memcpy(trGthHdr[numTr-1][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
			
			/*is this the highest xline found so far?*/
			if ((thdr.xline - xOrigin) > xNdxMx[thdr.iline - iOrigin]) xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;

			statusR = segyReadTrace(fin, &bh, &thdr, tempTr, ns, endian);
			
			if (thdr.iline > iMax){
				printf("\n\nFound iline greater than iMax.\n\n");
				return 0;
			}else if (thdr.iline < iOrigin){
				printf("\n\nFound iline less than iOrigin.\n\n");
				return 0;
			}else if (thdr.xline - xOrigin < 0){
				printf("\n\nFound xline  less than xOrigin.\n\n");
				return 0;
			}else if (thdr.xline > xMax){
				printf("\n\nFound xline greater than xMax.\n\n");
				return 0;
			}

		}
		
		/*I've stopped reading an inline.  If I have not stopped due to bad read (statusR !=0), then
		 * I have stopped due to finding the next inline.  Set the xNdxMn for this new iline*/
		if (statusR == 0){
			xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			ilineOld = thdr.iline;
		}
		/*read in new inline, time to update ilineCnt, inline counter/tracker*/
		ilineCnt++;
		//printf("\n\n ****************");
		//printf("\n inline %d", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline);

		sembl(trGth, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,thresh,xMax, rate, reso);
		
		if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
			printf("\n\n ****************");
			time(&now);
			printf("\n Finished Inline: %d Time: %s Percent Done: %d", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline - iOrigin + 1)/num10Percent*10 );
		
		}
		
		/*Output data - uncomment other for loop to output a padded version*/
		for (i = xNdxMn[ilineCnt]; i<=xNdxMx[ilineCnt];++i){ 
		//for (i = 0; i <= xMax - xOrigin;++i){ 	
			/*if no header exists (example padded trace within inline with iline = xline =0) put inline and xline numbers into headers*/
			if((*trGthHdr[numTr/2][i]).iline == 0){
				(*trGthHdr[numTr/2][i]).iline=(*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[numTr/2][i]).xline=xOrigin+i;
			}			
			segyWriteTrace(foutI, trGthHdr[numTr/2][i], dipTrI[numTr/2][i], ns, endian);
			segyWriteTrace(foutX, trGthHdr[numTr/2][i], dipTrX[numTr/2][i], ns, endian);
			segyWriteTrace(fsemb, trGthHdr[numTr/2][i], semTr[numTr/2][i], ns, endian);
		}	
	}

	/*ive exited (most likey due to end of file) while loop - print the last numTr/2 traces*/
	for (i=0;i<numTr/2;++i){
		for(j=xNdxMn[ilineCnt+i + 1];j<=xNdxMx[ilineCnt + i + 1];++j){
			
		//for(j=0;j<=xMax-xOrigin;++j){
			/*if no header exists (example padded trace within inline with iline = xline =0) put inline and xline numbers into headers*/
			if((*trGthHdr[numTr/2 + 1 + i][j]).iline == 0){
				(*trGthHdr[numTr/2 + 1 + i][j]).iline=(*trGthHdr[numTr/2 + 1 + i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[numTr/2 + 1 + i][j]).xline=xOrigin+j;
			}
			segyWriteTrace(foutI, trGthHdr[numTr/2 + 1 +i][j], dipTrI[numTr/2][j], ns, endian);
			segyWriteTrace(foutX, trGthHdr[numTr/2 + 1 + i][j], dipTrX[numTr/2][j], ns, endian);
			segyWriteTrace(fsemb, trGthHdr[numTr/2+ 1 + i][j], semTr[numTr/2 ][j], ns, endian);
		}
	}
	
printf("\n\nDONE!\n\n");
free(semTr);
free(trGth);
free(trGthHdr);
free(dipTrI);
free(dipTrX);
free(tempTr);
free(xNdxMx);
free(xNdxMn);
return 0;
}
   
/*
************************************************************************
*/
void doMessage(char *str)
{
   fprintf(stderr, "%s\n", str);
}
/*
***********************************************************************
*/

/***********************************************************************
 * Function: Windowed Semblance along a plane in 3D
 * Input: Gather of traces, current inline and xline of analysis, 
 * 		  current sample point, sembalance window size, and max dip
 * 
 * Description:  Performs windowed semblance, where the semblance is
 * 				 computed as the power of the resulting stacked trace
 *  			 divided by the sum of the power of the indvidual traces
***********************************************************************/
void sembl(float ***trGth,int ns, int numTr,int iNdx ,int window,int GlobalmaxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx, float thresh, int xMax, int rate, int reso){

int thread_id, nloops;

	float **squareStkA; //array of squared stk values for each value of dipX and dipY centered about sample
	float **squares; //array of squared sample values for each value of dipX and dipY centered about sample

#pragma omp parallel private(thread_id, nloops, squareStkA, squares) 
{
	int i,j,k;
	
	int maxDip;
	maxDip=GlobalmaxDip*reso;
	
	float stack=0.0; //holds the stacked values for each value of window
	float squareStk=0; //sum of squares and sum of stcked squares
	float sembVal=0, hold=0, holdDip=0;
	int x=0,y=0; //
	int slopeX; //slope in xline direction 
	int slopeY; //slope in inline direction
	int dipIl, dipXl; //dip for xline and inline direction
	int samp = 0;
	int xStart = xNdxMn[iNdx];//-xOrigin;//index number of first xline.  if xNdxMn = xOrigin, then xStart =0 
	int xEnd = xNdxMx[iNdx];//-xOrigin; //xEnd indx number of last xline
	int xNdx;
	int shift = (window/2 + 2*numTr/2*maxDip +10); //must shift down for semblace calculation 	
	float subSquares = 0; //amnt to subract from squares for that particular dipX and dipY
	
	float subStk = 0.0; // amnt to subtract from sqaureStkA for that particular dipX and dipY
	int halfTr = numTr / 2;
	int halfWin = window / 2;
	int totalSlope = 0;
	int totalSlopePlus =0;
	float traceCalc;
	float subTraceCalc;
	float scale;
	float holdDipX=0;
	float holdDipY=0;
	int boundXlow, boundXHigh;
		
	float **squareStkA; //array of squared stk values for each value of dipX and dipY centered about sample
	float **squares; //array of squared sample values for each value of dipX and dipY centered about sample
	
	/*if the resolution is 2 per sample and you want 2 samples then you gotta calculate 2*2 */
	
	/*Allocating some memory*/
	squareStkA = (float**)calloc((2*maxDip+8) , sizeof(float *));
    for (i=0; i<(2*maxDip+8); ++i){
         squareStkA[i] = (float *)calloc((2*maxDip+8),sizeof(float));   
	 }
	squares = (float**)calloc((2*maxDip+8) , sizeof(float *));
    for (i=0; i<(2*maxDip+8); ++i){
         squares[i] = (float *)calloc((2*maxDip+8),sizeof(float));    
	}
	
/* iNdx is the array index number for the inline direction of trGth*/	
/* xNdx is the array index number for the xline direction*/
	
	/*walking down the inline.  Remeber, this will start index numTr/2*/
	#pragma omp for
	//for (xNdx=xStart+(numTr/2);xNdx <= xEnd-(numTr/2); ++xNdx){
	for (xNdx=xStart;xNdx <= xEnd; ++xNdx){
		
		/*setting bounds*/
		if(xNdx-xStart < halfTr){
			boundXlow=xStart - xNdx + halfTr;
		}else {
			boundXlow= 0;
		}
		if( xEnd - xNdx < halfTr){
			boundXHigh= xNdx -xEnd + halfTr;
		}else {
			boundXHigh = 0;
		}


			/*computing the first sample point, then and window values, then i will loop through the rest*/
			samp=shift;				
			/*loop thru all dips 1 sample at a time in both inline and xline direction*/
			for (dipIl= -maxDip;dipIl<=maxDip;dipIl++){
				for(dipXl = -maxDip;dipXl<=maxDip;++dipXl){

					/*reset arrays for next xline trace*/
					squareStkA[maxDip + dipIl][maxDip + dipXl] = 0.0;
					squares[maxDip + dipIl][maxDip + dipXl] = 0.0;
					
						/*loop though the window and gather to produce a slant stack centered
						* about the center trace of gather and window*/	
						for (k=0; k<window; ++k){
							/*for the number of traces*/
							for(x = -(halfTr ); x <= (halfTr); ++x){
								
								for (y =  -(halfTr - boundXlow); y <= (halfTr - boundXHigh); ++y){

									/*slope calcultion starts in upper left corner. Remember, dip 
									* increments starting as a NEGATIVE number*/
									slopeX=(dipIl*(x));
									slopeY = (dipXl*(y));

									scale = fabs(((float)((abs(-slopeX - slopeY))%rate))/rate);

									totalSlope = samp + (- slopeX - slopeY)/reso;
							
									if(slopeX+slopeY > 0){
										totalSlopePlus = (float)totalSlope - 1.05;
									}else if(slopeX + slopeY < 0){
										totalSlopePlus = (float)totalSlope + 1.05;
									}else {
										totalSlopePlus = (float)totalSlope;
									}
									//#pragma omp critical 
									traceCalc= (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin +k] + (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin +k];
									
									stack+=traceCalc;
									squares[maxDip + dipIl][maxDip + dipXl] += traceCalc*traceCalc;

									/*holding onto first value to subtract from values later*/
									if(k==0){
										subStk+=traceCalc;
									}
								}
							}
							
					//	}// End of loop
						
							/*subSquares value to be removed later*/
							if (k == 0){
								subSquares = squares[maxDip + dipIl][maxDip + dipXl];
							}
					
						/*gettting the square of the stack and then zeroing stack value for later use*/
						squareStkA[maxDip+dipIl][maxDip+dipXl]+=stack*stack;
						//printf("\n\n test 2");
						stack=0.0;
						}
						/*computing semblance value*/
						if(squares[maxDip+dipIl][maxDip+dipXl] < .000000001){
							sembVal =0;
						}else{
						sembVal=squareStkA[maxDip+dipIl][maxDip+dipXl]/((squares[maxDip+dipIl][maxDip+dipXl])*numTr*numTr);
						}

						/*is this the best semblance value yet? and is it at least as good as thresh*/
						
						if(sembVal > hold){
							if(sembVal >= thresh){
								dipTrI[halfTr][xNdx][samp] = dipIl*rate/reso;
								dipTrX[halfTr][xNdx][samp] = dipXl*rate/reso;
								holdDipX = dipIl;
								holdDipY = dipXl;
							}else{
								dipTrI[halfTr][xNdx][samp] = 0;
								dipTrX[halfTr][xNdx][samp] = 0;
							}
							semTr[halfTr][xNdx][samp] = 1.0 - sembVal;
							hold = sembVal;
						}
				/*by now I have computed the semblance for this particular dip. Before I go further
				I will subtract the values assiciated with the value k=0 (the first values in the window) from this dip*/
				squareStkA[maxDip+dipIl][maxDip+dipXl] -= subStk*subStk;
				squares[maxDip+dipIl][maxDip+dipXl]-= subSquares;
				subSquares=0.0;
				subStk=0.0;
				}	
			}
			hold = 0.0;
			holdDipX=holdDipY=0;

		/*setting k = window-1,I only need to add the value for the last sample in the window
		* back into squareStkA for semblance calculation*/
		k = window - 1;
		
		/*running through remaining samples*/
		for(samp=shift+ 1;samp<=ns-shift; ++samp){	  	
			for (dipIl= -maxDip;dipIl<=maxDip;++dipIl){
				for(dipXl = -maxDip;dipXl<=maxDip;++dipXl){

					/*running through the number of traces*/
						for (x = -(halfTr); x <= (halfTr); ++x){
							//for (y = -(halfTr); y <= (halfTr); ++y){
							 for (y =  -(halfTr - boundXlow); y <= (halfTr - boundXHigh); ++y){
								 
								/*slope calcultion starts in upper left corner. Remember, dip 
								* increments starting as a NEGATIVE number*/
								slopeX=(int)(dipIl*(x) );
								slopeY=(int)(dipXl*(y) );

								scale = fabs(((float)((abs(-slopeX - slopeY))%rate))/rate);

								totalSlope = samp + (- slopeX - slopeY)/reso;
							
									if(slopeX+slopeY > 0){
										totalSlopePlus = (float)totalSlope - 1.05;
									}else if(slopeX + slopeY < 0){
										totalSlopePlus = (float)totalSlope + 1.05;
									}else {
										totalSlopePlus = (float)totalSlope;
									}
							//	#pragma omp critical 
								traceCalc = (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin +k] + (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin +k];
							//	#pragma omp critical 
								subTraceCalc = (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin ]+ (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin ];
								
								stack += traceCalc;
								subStk += subTraceCalc; //subtract out the k from totalSlope to compute value of window for k=0
								squares[maxDip + dipIl][maxDip + dipXl] += traceCalc*traceCalc;
								subSquares += subTraceCalc*subTraceCalc;	
							}
						}
						
						squareStkA[maxDip+dipIl][maxDip+dipXl]+= stack*stack;
						stack = 0.0;
						
						if(squares[maxDip+dipIl][maxDip+dipXl] < .000000001){
							sembVal =0;
						}else{
						sembVal=squareStkA[maxDip+dipIl][maxDip+dipXl]/((squares[maxDip+dipIl][maxDip+dipXl])*numTr*numTr);
						}
						
						if(sembVal >=hold){
							if(sembVal >= thresh){
								dipTrI[halfTr][xNdx][samp] = dipIl*rate/reso;
								dipTrX[halfTr][xNdx][samp] = dipXl*rate/reso;
								holdDipX = dipIl;
								holdDipY = dipXl;
								
							}else{
								dipTrI[halfTr][xNdx][samp] = 0;
								dipTrX[halfTr][xNdx][samp] = 0;
							}
							semTr[halfTr][xNdx][samp] = 1.0 - sembVal;
							hold = sembVal;
						}
						/*once the semblance value has been calculated, I can subtract out subSquares and subStk*/
						squareStkA[maxDip + dipIl][maxDip + dipXl] -= subStk*subStk;
						subStk = 0.0;
						squares[maxDip + dipIl][maxDip + dipXl] -= subSquares;
						subSquares = 0.0;
					}
				}
			hold = 0.0;
			holdDipX=holdDipY=0;
		}
		
	}//end of xndx
free(squareStkA);
free(squares);
}//endo of openmp

}


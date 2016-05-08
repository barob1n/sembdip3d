#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defns.h"
#include "segy.h"
#include "segyIO_class.h"
#include <time.h>

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
 
void sembl(float ***trGth, int ns, int numTr, int iNdx, int window, int maxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, float thresh, int xMax, int rate);

int main(int argc, char *argv[])
{
FILE  *fin=NULL,  *foutX=NULL, *foutI=NULL, *fsemb=NULL; //fout = file out , fmodel = model input
int   fmt, i, itr, ierr, k;
int   ns;	/* Number of samples in a trace	*/
float dts;	/* Sample rate in seconds	*/
segy  thdr;	/* Segy trace header		*/
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

float temp;
endian = checkEndian();

/* The command line format is avoProg <output.sgy> <model.txt>*/
	if(argc ==13){
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
		if(sscanf(argv[9], "%d", &iOrigin) != 1){
			printf("\niOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[10], "%d", &xOrigin) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[11], "%d", &iMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[12], "%d", &xMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}		
	} else {
	fprintf(stderr,"\n***************************************************************************\n");
      fprintf(stderr,"Computes the dip and coherancy at each sample via semblance. \n\n");
      fprintf(stderr,"Program expects the following command line: \n ");
      fprintf(stderr,"sembDip <in.sgy><dipI.sgy><dipX.sgy><sem.sgy><Num><maxDip><win><thresh>\n");
      fprintf(stderr,"<iOrigin><xOrigin><iMax><xMax>\n\n");
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
      fprintf(stderr,"iOrigin: Starting inline number.\n\n");
      fprintf(stderr,"xOrigin: Starting xnline number.\n\n");
      fprintf(stderr,"iMax: Max iline. \n\n");
      fprintf(stderr,"xMax: Max xline.");
    fprintf(stderr, "\n**************************************************************************\n");
      return 0;
    } 

	/*print parameters. ask user to verify they are correct-option to abort*/
	printf("\nInput paramters are:");
	printf("\n\n  Input seismic: %s", argv[1]);
	printf("\n  Output dipI: %s", argv[2]);
	printf("\n  Output dipX: %s", argv[3]);
	printf("\n  Output semb: %s", argv[4]);
	//printf("\n Output velocity:  %s", fvelout);
	printf("\n  Number to trace for semb: %d", numTr);
	printf("\n  Max dip (in samples!): %d", maxDip);
	printf("\n  Window (in samples!): %d", window);
	printf("\n  threshhold: %f", thresh);
	printf("\n  inline origin: %d", iOrigin);
	printf("\n  xline origin: %d", xOrigin);
	printf("\n  iMax : %d", iMax);
	printf("\n  xMax: %d", xMax);
	
	printf("\n\nAre these correct? (-1: exit)");
	fgets(buff,80,stdin);
	sscanf(buff,"%d",&choice);
	if(choice == -1){ return 0;} 
	
	segyReadHeader(fin, chdr, &bh, endian);
	ns=bh.hns;
	rate = bh.hdt;
	rate=rate/1000;
	printf("\nNumber of samples: %d", ns);
	//printf("\n format: %d\n\n", bh.format);
	
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
	/*are the values sensabale?*/
	if (thdr.iline > iMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}else if (thdr.xline - xOrigin < 0){
		printf("\n\nFirst xline read less than xOrigin.\n\n");
		return 0;
	}else if (thdr.xline > xMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}
	
	/*save the xline values and set iline switch -really only need to set the Mn value*/
	xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	ilineOld=thdr.iline;
	
	/*read in teh first numTr number of traces*/
	while ( thdr.iline - iOrigin < numTr){
	
		/*copy the values read from previous trace read*/
		memcpy(trGth[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTr,(ns+1)*sizeof(float));
		memcpy(trGthHdr[thdr.iline - iOrigin][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
		
		statusR = segyReadTrace(fin, &bh, &thdr,tempTr, ns, endian);	
			
		/*if it gets to end of reading before having read in enough traces abort
		*else end o file and you have enough lines exit while loop */
		if(statusR!=0){
			if(thdr.iline - iOrigin < numTr -1 ){
				printf("\nNot enough ilines for a single semblance calc!\n");
				return  0;
			}else{
					
				break; /*minimum number of lines met...exit loop*/
			}
		}
		
		    /*the xline number just read is large than xNdxMx*/
			if ( (thdr.xline-xOrigin) > xNdxMx[thdr.iline-iOrigin]) xNdxMx[thdr.iline - iOrigin]=thdr.xline - xOrigin;
			if (ilineOld != thdr.iline){
				xNdxMn[thdr.iline-iOrigin]=thdr.xline-xOrigin;
				ilineOld=thdr.iline;
				//printf("\n\n I found %d\n\n", thdr.iline);
			}
	}

	/*copy the value read last, which caused while loop to exit*/
	ilineCnt=numTr/2;
	printf("\n\n ****************");
	printf("\n inline %d", thdr.iline - numTr/2 -1);

	/*compute semblance and dips*/
	sembl(trGth, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,xMax,rate);

	/*spits out first numTr/2 +1  lines*/
	for (i=0;i<=numTr/2;++i){
		//for(j=xNdxMn[ilineCnt-numTr/2 + i];j<=xNdxMx[ilineCnt-numTr/2 + i];++j){
		for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
				segyWriteTrace(foutI, trGthHdr[i][j], dipTrI[numTr/2][j], ns, endian);
				segyWriteTrace(foutX, trGthHdr[i][j], dipTrX[numTr/2][j], ns, endian);
				segyWriteTrace(fsemb, trGthHdr[i][j], semTr[numTr/2][j], ns, endian);
			}
			//printf("\n\n trace header first :%d\n\n", (*trGthHdr[i][0]).iline);
	}
	
   while (statusR == 0){ 

		/*Shift all the traces back by one, meaning that
		* trGth[1] -> trGth[0], trGth[2] -> [1], and tempTr -> trGth[numTr-1]*/
		for (k =0 ;k<numTr-1;++k){
			for(j=0;j<xMax-xOrigin+1;++j){
					memcpy(trGth[k][j],trGth[k+1][j],(ns+1)*sizeof(float));
					memcpy(trGthHdr[k][j],trGthHdr[k+1][j], sizeof(*trGthHdr[0][0]));
			}
		}
		/*getting another inline*/
		/*setting last inline to zero before I read a new inline into it*/
		for (j =0; j < xMax-xOrigin +1;++j){
			memset(trGth[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(trGthHdr[numTr-1][j],0,sizeof(*trGthHdr[0][0]));
		}

		/*Copying the last read tempTr and and thdr before continuing with rest of line.  This is done 
		since the program only know when to stop reading an inline once it has read in the first xline 
		of the next inline.  So, this value is being added back into the gathger*/
	
		//memcpy(trGth[numTr-1][thdr.xline -xOrigin], tempTr,(ns+1)*sizeof(float));
		//memcpy(trGthHdr[numTr-1][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
		while(statusR ==0 && thdr.iline ==ilineOld){
			
			memcpy(trGth[numTr - 1][thdr.xline - xOrigin], tempTr, (ns + 1)*sizeof(float));
			memcpy(trGthHdr[numTr-1][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));

			if ((thdr.xline - xOrigin) > xNdxMx[thdr.iline - iOrigin]) xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;

			statusR = segyReadTrace(fin, &bh, &thdr, tempTr, ns, endian);
		}

		if (statusR == 0){
			xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			ilineOld = thdr.iline;
		}

		ilineCnt++;//starts at second inline, ;
		printf("\n\n ****************");
		printf("\n inline %d", thdr.iline - numTr/2 -1);

		sembl(trGth, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,xMax, rate);
		
		 /*output inline*/
		for (i = xNdxMn[ilineCnt]; i<=xNdxMx[ilineCnt];++i){ 			
			segyWriteTrace(foutI, trGthHdr[numTr/2][i], dipTrI[numTr/2][i], ns, endian);
			segyWriteTrace(foutX, trGthHdr[numTr/2][i], dipTrX[numTr/2][i], ns, endian);
			segyWriteTrace(fsemb, trGthHdr[numTr/2][i], semTr[numTr/2][i], ns, endian);
		}
	}

	/*ive exited (most likey due to end of file) while loop - print the last numTr/2 traces*/
	for (i=0;i<numTr/2;++i){
		//printf("\n\n TEST last lines \n\n ");
		for(j=xNdxMn[ilineCnt+i];j<=xNdxMx[ilineCnt+i];++j){
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
/*void medFilt(float ***dipTrI, float ***dipTrX, float ***semTr,int *xNdxMn, int *xNdxMx int numTr,int ns){
	
	int i,j,k,samp;
	int shift = (window/2 + 2*(numTr/2)*maxDip);
	
	
	for(i=xNdxMn+1;i<xNdxMx){
		for (samp =shift; samp <= ns-shift;++samp){
			
}*/
/***********************************************************************
 * Function: Windowed Semblance along a plane in 3D
 * Input: Gather of traces, current inline and xline of analysis, 
 * 		  current sample point, sembalance window size, and max dip
 * 
 * Description:  Performs windowed semblance, where the semblance is
 * 				 computed as the power of the resulting stacked trace
 *  			 divided by the sum of the power of the indvidual traces
***********************************************************************/
void sembl(float ***trGth,int ns, int numTr,int iNdx ,int window,int maxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, float thresh, int xMax, int rate){
	
	int i,j,k;
	float stack=0.0; //holds the stacked values for each value of window
	float squareStk=0; //sum of squares and sum of stcked squares
	float sembVal=0, hold=0, holdDip=0;
	int x=0,y=0; //
	int slopeX; //slope in xline direction 
	int slopeY; //slope in inline direction
	int dzX, dzY; //dip for xline and inline direction
	int samp  =0;
	int xStart=xNdxMn[iNdx];//-xOrigin;//index number of first xline.  if xNdxMn = xOrigin, then xStart =0 
	int xEnd = xNdxMx[iNdx];//-xOrigin; //xEnd indx number of last xline
	int xNdx;
	int shift = (window/2 + 2*(numTr/2)*maxDip); //must shift down for semblace calculation 
	float subSquares=0; //amnt to subract from squares for that particular dipX and dipY
	float **squareStkA; //array of squared stk values for each value of dipX and dipY centered about sample
	float **squares; //array of squared sample values for each value of dipX and dipY centered about sample
	float subStk = 0.0; // amnt to subtract from sqaureStkA for that particular dipX and dipY
	int halfTr = numTr / 2;
	int halfWin = window / 2;
	int totalSlope = 0;
	int totalSlopePlus =0;
	float traceCalc;
	float subTraceCalc;
	float scale;
	int test;
	float holdDipX=0;
	float holdDipY=0;
	int a=5, b=5;
	time_t now;

	/*Allocating some memory*/
	squareStkA = (float**)calloc((2*maxDip+2) , sizeof(float *));
    for (i=0; i<(2*maxDip+2); ++i){
         squareStkA[i] = (float *)calloc((2*maxDip+2),sizeof(float));   
	 }
	squares = (float**)calloc((2*maxDip+2) , sizeof(float *));
    for (i=0; i<(2*maxDip+2); ++i){
         squares[i] = (float *)calloc((2*maxDip+2),sizeof(float));    
	}
/* iNdx is the array index number for the inline direction of trGth*/	
/* xNdx is the array index number for the xline direction*/
	time(&now);
	 printf("\n xStart %d  ", xStart); 
	 printf("\n xEnd %d  ", xEnd);
	 printf("\n Time: %s", ctime(&now));
	 
	/*walking down the inline.  Remeber, this will start index numTr/2*/
	for (xNdx=xStart+(numTr/2);xNdx <= xEnd-(numTr/2); ++xNdx){
		
			/*computing the first sample point, then and window values, then i will loop through the rest*/
			samp=shift;				
			/*loop thru all dips 1 sample at a time in both inline and xline direction*/
			for (dzX= -maxDip;dzX<=maxDip;dzX++){
				for(dzY = -maxDip;dzY<=maxDip;++dzY){

					/*reset arrays for next xline trace*/
					squareStkA[maxDip + dzX][maxDip + dzY] = 0.0;
					squares[maxDip + dzX][maxDip + dzY] = 0.0;
					
						/*loop though the window and gather to produce a slant stack centered
						* about the center trace of gather and window*/	
						for (k=0; k<window; ++k){

							/*for the number of traces*/
							for(x = -(halfTr); x <= (halfTr); ++x){
								
								for (y = (-halfTr); y <= (halfTr); ++y){

									/*slope calcultion starts in upper left corner. Remember, dip 
									* increments starting as a NEGATIVE number*/
									slopeX=(dzX*(x));
									slopeY = (dzY*(y));
									
									//scale = fabs(((float)((rate -slopeX - slopeY)%rate))/rate);
									scale = fabs(((float)((abs(-slopeX - slopeY))%rate))/rate);
									
									//totalSlope = (samp*rate - slopeX - slopeY)/rate;
									totalSlope = samp + (- slopeX - slopeY)/rate;
							
									if(slopeX+slopeY > 0){
										totalSlopePlus = (float)totalSlope - 1.05;
									}else if(slopeX + slopeY < 0){
										totalSlopePlus = (float)totalSlope + 1.05;
									}else {
										totalSlopePlus = (float)totalSlope;
									}
									traceCalc= (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin +k] + (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin +k];
									
									stack+=traceCalc;
									squares[maxDip + dzX][maxDip + dzY] += traceCalc*traceCalc;
									
								//printf("\n\n Current Sembval %f", stack*stack/(.0000001+squares[maxDip + dzX][maxDip + dzY]));
								//printf("\n\n big Lop448  %d, %f %f %f:", samp, scale,traceCalc, subTraceCalc);
								//printf("\n\n %f % f % f % f ", stack, subStk, squares[maxDip + dzX][maxDip + dzY], subSquares);
								//getchar();
									
									
									/*holding onto first value to subtract from values later*/
									if(k==0){
										subStk+=traceCalc;
									}
								}
							}
							/*subSquares value to be removed later*/
							if (k == 0){
								subSquares = squares[maxDip + dzX][maxDip + dzY];
							}
						
						/*gettting the square of the stack and then zeroing stack value for later use*/
						squareStkA[maxDip+dzX][maxDip+dzY]+=stack*stack;
						stack=0.0;
						}
						/*computing semblance value*/
						if(squares[maxDip+dzX][maxDip+dzY] < .000000001){
							sembVal =0;
						}else{
						sembVal=squareStkA[maxDip+dzX][maxDip+dzY]/((squares[maxDip+dzX][maxDip+dzY])*numTr*numTr);
						}
						//printf("\n\n sembVale %f:", sembVal);
						//getchar();

						/*is this the best semblance value yet? and is it at least as good as thresh*/
						if(sembVal > hold){
							if(sembVal >= thresh){
								dipTrI[halfTr][xNdx][samp] = dzX;
								dipTrX[halfTr][xNdx][samp] = dzY;
								holdDipX = dzX;
								holdDipY = dzY;
							}else{
								dipTrI[halfTr][xNdx][samp] = 0;
								dipTrX[halfTr][xNdx][samp] = 0;
							}
							semTr[halfTr][xNdx][samp] = 1.0 - sembVal;
							hold = sembVal;
						}
				/*by now I have computed the semblance for this particular dip. Before I go further
				I will subtract the values assiciated with the value k=0 (the first values in the window) from this dip*/
				squareStkA[maxDip+dzX][maxDip+dzY] -= subStk*subStk;
				squares[maxDip+dzX][maxDip+dzY]-= subSquares;
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
			for (dzX= -maxDip;dzX<=maxDip;++dzX){
				for(dzY = -maxDip;dzY<=maxDip;++dzY){
					
					//squares[maxDip + dzX][maxDip + dzY] = 0.0;
			
					/*running through the number of traces*/
						for (x = -(halfTr); x <= (halfTr); ++x){
							for (y = -(halfTr); y <= (halfTr); ++y){

								/*slope calcultion starts in upper left corner. Remember, dip 
								* increments starting as a NEGATIVE number*/
								slopeX=(int)(dzX*(x) );
								slopeY=(int)(dzY*(y) );

								//scale = fabs(((float)((rate -slopeX - slopeY)%rate))/rate);
								scale = fabs(((float)((abs(-slopeX - slopeY))%rate))/rate);
								
								//totalSlope = (float)samp - (float)slopeX - (float)slopeY;
								totalSlope = samp + (- slopeX - slopeY)/rate;
							
									if(slopeX+slopeY > 0){
										totalSlopePlus = (float)totalSlope - 1.05;
									}else if(slopeX + slopeY < 0){
										totalSlopePlus = (float)totalSlope + 1.05;
									}else {
										totalSlopePlus = (float)totalSlope;
									}
								traceCalc = (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin +k] + (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin +k];
								
							
								subTraceCalc = (1.0-scale)*trGth[halfTr+x][xNdx+y][totalSlope - halfWin ]+ (scale)*trGth[halfTr+x][xNdx+y][totalSlopePlus - halfWin ];
								
								stack += traceCalc;
								subStk += subTraceCalc; //subtract out the k from totalSlope to compute value of window for k=0
								squares[maxDip + dzX][maxDip + dzY] += traceCalc*traceCalc;
								subSquares += subTraceCalc*subTraceCalc;
								
							}
						}
						
						squareStkA[maxDip+dzX][maxDip+dzY]+= stack*stack;
						stack = 0.0;
						if(squares[maxDip+dzX][maxDip+dzY] < .000000001){
							sembVal =0;
						}else{
						sembVal=squareStkA[maxDip+dzX][maxDip+dzY]/((squares[maxDip+dzX][maxDip+dzY])*numTr*numTr);
						}
						
						if(sembVal >=hold){
							if(sembVal >= thresh){
								dipTrI[halfTr][xNdx][samp] = dzX;
								dipTrX[halfTr][xNdx][samp] = dzY;
								holdDipX = dzX;
								holdDipY = dzY;
								
							}else{
								dipTrI[halfTr][xNdx][samp] = 0;
								dipTrX[halfTr][xNdx][samp] = 0;
							}
							semTr[halfTr][xNdx][samp] = 1.0 - sembVal;
							hold = sembVal;
						}
						/*once the semblance value has been calculated, I can subtract out subSquares and subStk*/
						//printf("\n sembValue final  %f:", hold);
						//getchar();
						squareStkA[maxDip + dzX][maxDip + dzY] -= subStk*subStk;
						subStk = 0.0;
						squares[maxDip + dzX][maxDip + dzY] -= subSquares;
						subSquares = 0.0;
					}
				}
			hold = 0.0;
			holdDipX=holdDipY=0;
		}
}
free(squareStkA);
free(squares);
}
/***********************************************************************
 * Function: Create Array 
 * Input: array name, x,y,z dimensions.
 * 
 * Description:  creates pointer array.
***********************************************************************/


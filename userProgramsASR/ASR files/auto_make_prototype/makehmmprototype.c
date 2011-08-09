/* MakeHMMPrototype.c */

#include <stdio.h>
#include <stdlib.h>

void printUsage(char * arg0) {
  printf("%s, written by David Dean, 2004\n\n",arg0);
  printf("Generates a forward only, diagonal covariance HMM prototype\n\n");
  printf("usage:\n");
  printf("   %s name parameterkind vectorsize numberstates mixturesperstate \\\n",arg0);
  printf("                [ streamcount vecsize1 ... vecsizeN [ weight1 ... weightN ] ]\n\n");
  printf("   name             : name of hmm (eg. prototype)\n");
  printf("   parameterkind    : parameter kind (eg. MFCC_0_D_A)\n");
  printf("   vectorsize       : size of feature vector\n");
  printf("   numberstates     : total number of states\n");
  printf("   mixturesperstate : number of mixtures per state\n");
  printf("   streamcount ...  : multiple-stream support\n");
  printf("      vecsizeN      : size of each stream (must total to vectorsize)\n");
  printf("      weightN       : weight of each stream on final score (defaults to 1)\n");
}  

/* Main Function */
int main (int argc, char *argv[]) {
  if (argc < 6) { /* 5+ args + commandline */
    printUsage(argv[0]);
    exit(1);
  }
  
  char* NAME = argv[1];
  char* PARAMETER_KIND = argv[2];
  int VECTOR_SIZE = atoi(argv[3]);
  int NUMBER_STATES = atoi(argv[4]);
  int MIXTURES_PER_STATE = atoi(argv[5]);
  
  int NUMBER_STREAMS = 1;
  int * STREAM_VECTOR_SIZE = 0;
  double * STREAM_WEIGHT = 0;
  
  float MIXTURE_WEIGHT = 1.0 / MIXTURES_PER_STATE;
  float INITIAL_MEAN = 0.0;
  float INITIAL_VARIANCE = 1.0;
  
  int stream, state, mixture, feature, statefrom, stateto;
  
  /* load stream info from command line */
  
  if (argc > 6) { /* number of streams specified on command line */
    NUMBER_STREAMS = atoi(argv[6]);
    STREAM_VECTOR_SIZE = (int *) malloc( sizeof(int) * NUMBER_STREAMS );
    STREAM_WEIGHT = (double *) malloc( sizeof(double) * NUMBER_STREAMS );
    for (stream = 0; stream < NUMBER_STREAMS; stream++) {
      STREAM_VECTOR_SIZE[stream] = atoi(argv[7+stream]);
      if (argc > 7 + NUMBER_STREAMS) { /* weights of streams specified on cmd line */ 
	STREAM_WEIGHT[stream] = atof(argv[7+NUMBER_STREAMS+stream]);
      }	else {
	STREAM_WEIGHT[stream] = 1.0;
      }
    }
  } else { /* one stream of full vector size */
    NUMBER_STREAMS = 1;
    STREAM_VECTOR_SIZE = (int *) malloc(sizeof(int) * 1);
    STREAM_WEIGHT = (double *) malloc( sizeof(float) * 1);
    STREAM_VECTOR_SIZE[0] = VECTOR_SIZE;
    STREAM_WEIGHT[0] = 1.0;
  }    
  
  /* output header */
  printf("~o <VecSize> %d <%s>\n",VECTOR_SIZE,PARAMETER_KIND);
  if (NUMBER_STREAMS > 1) {
    printf("   <StreamInfo> %d",NUMBER_STREAMS);
    for (stream = 0; stream < NUMBER_STREAMS; stream++) {
      printf(" %d",STREAM_VECTOR_SIZE[stream]);
    }
    printf("\n");
  }
  printf("~h \"%s\"\n",NAME);
  printf("<BeginHMM>\n");
  printf(" <NumStates> %d\n",NUMBER_STATES);
  
  /* output each emitting state (ie. not the first or last) */
  for (state=2; state<NUMBER_STATES; state++) {
    
    printf(" <State> %d\n",state);
    
    printf("  <NumMixes>"); 
    for (stream=1; stream<=NUMBER_STREAMS; stream++) {
      printf(" %d",MIXTURES_PER_STATE);
    }
    printf("\n");
     
    printf("  <SWeights> %d",NUMBER_STREAMS);
    for (stream=1; stream<=NUMBER_STREAMS; stream++) {
      printf(" %g",STREAM_WEIGHT[stream-1]);
    }
    printf("\n");

    for (stream=1; stream<=NUMBER_STREAMS; stream++) {

      if (NUMBER_STREAMS > 1) {
	printf("  <Stream> %d\n",stream);
      }

      /* output each mixture for each state */
      for (mixture=1; mixture<=MIXTURES_PER_STATE; mixture++) {
	printf("   <Mixture> %d %g\n",mixture,MIXTURE_WEIGHT);
	
	/* output mean values for each mixture */
	printf("    <Mean> %d\n",STREAM_VECTOR_SIZE[stream-1]);
	printf("    ");
	for (feature=1; feature<=STREAM_VECTOR_SIZE[stream-1]; feature++)
	  {
	    printf(" %f",INITIAL_MEAN);
	  }
	printf("\n");
	
	/* output diagonal covariance values for each mixture */
	printf("    <Variance> %d\n",STREAM_VECTOR_SIZE[stream-1]);
	printf("    ");
	for (feature=1; feature<=STREAM_VECTOR_SIZE[stream-1]; feature++) {
	  printf(" %f",INITIAL_VARIANCE);
	}
	printf("\n");
      }
    }
  }
  
  /* transition matrix */
  printf(" <TransP> %d\n",NUMBER_STATES);
  for (statefrom=1; statefrom<=NUMBER_STATES; statefrom++)
    {
      printf(" ");
      for (stateto=1; stateto<=NUMBER_STATES; stateto++)
	{
	  if (statefrom==1 && stateto==2) 
	    printf(" %f",1.0);
	  else if (statefrom>1 && statefrom<NUMBER_STATES && (stateto==statefrom || stateto==(statefrom+1))) 
	    printf(" %f",0.5);
	  else 
	    printf(" %f",0.0);
	}
      printf("\n");
    }
  
  /* output footer */
  printf("<EndHMM>\n");

  /* free memory from mallocs */
  if ( STREAM_VECTOR_SIZE ) free( STREAM_VECTOR_SIZE );
  if ( STREAM_WEIGHT ) free( STREAM_WEIGHT );
}

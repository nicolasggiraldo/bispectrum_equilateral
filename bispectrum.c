#include       <stdio.h>
#include      <stdlib.h>
#include        <math.h>
#include       <fftw3.h>
#include      <string.h>
#include <gsl/gsl_rng.h>
#include        <time.h>
#include      <unistd.h>
#include      <limits.h>
#include         <mpi.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

/*
  FFTW_MEASURE: find the optimal plan by actually computing several FFTs  
  FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
  FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
  FFTW_IN_PLACE: a plan assumes that the in and out are same
*/

int main(int argc, char *argv[]){
  int i, j, k, l;                  // Counter in the X, Y, and Z axis.
  unsigned long int rand_i;        /* Index for the subsample of density contrast */
  unsigned long int rand_j;        /* in the range between                        */
  //unsigned long int rand_k; /* k1-DELTA_K/2 to k1+DELTA_K/2                */
  long int id_cell;                // Counter for the id of the cell
  double kMag;                     /* Magnitude of the wavenumber vector */
  FILE *fout=NULL;                 // File handler to output
  double *kpos;                    /* Array of positions values according to FFTW 
				      k-position convention */ 
  int *indexpos=NULL;              /* Array of index values according to FFTW 
				      k-position convention */
  int indexcord[2];
  int m3[3];
  int rank, size;
  int Nbins;
  int *taskBin=NULL;
  struct densityContrast *q1=NULL; /* Array of structure with the subset of 
				      density contrast with range between 
				      k1-DELTA_K/2 to k1+DELTA_K/2 */
  struct binStruct *bindata=NULL;  /* Data structure with information of the 
				      measures in the taken bins */
  
  MPI_Status status;
  int chunk;

  /*
  gsl_rng *r=NULL;
  long seed;
  const gsl_rng_type *T;
  seed = time(NULL) * getpid();
  T = gsl_rng_gfsr4;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  //*/
  
  
  
  /////////////////////
  //* MPI BEGGINING *//
  /////////////////////
  MPI_Init(&argc, &argv);
  
  // Number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Task numbering
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(rank == 0){
    if(argc<2){
      printf("\n***********************************");
      printf("***********************************\n");
      printf("%s: You must specify the name of the parameter file\n",argv[0]);
      printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  

  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1], rank) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad path (or name) to the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad settings in the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////
  //* READING CELL BINARY FILE *//
  ////////////////////////////////
  switch( readBinaryFile(rank) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: The parameter file could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: FFTW arrays could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////
  //* EXECUTING FFTW3 ROUTINES *//
  ////////////////////////////////
  if(rank == 0){
    /* Do forward FFT */
    fftw_execute(forwardPlan);
    printf("\n-----------------------------------------------\n");
    printf("Fourier Transform succes.\n");
    fftw_free(denConX);
    fftw_destroy_plan(forwardPlan);
  }

  //fftw_free(denConX);
  //fftw_destroy_plan(forwardPlan);
  
  /* Waiting the processor 0 to do the FFT */
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* Sending all the deltak grid to the others processors */
  chunk = GV.NGRID3/sizeof(fftw_complex);
  for(l=0; l<sizeof(fftw_complex); l++)
    MPI_Bcast(denConK+(l*chunk), chunk*sizeof(fftw_complex), MPI_BYTE, 0, MPI_COMM_WORLD); 
  
  
  
  ///////////////////////////////////////
  //* SETTING FFTW3 COORDINATE SYSTEM *//
  ///////////////////////////////////////
  
  /* Position array for storing in the densityContrast */
  kpos = (double *) calloc(GV.NGRID, sizeof(double));
  if(kpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: kpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  indexpos = (int *) calloc(GV.NGRID, sizeof(int));
  if(indexpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: indexpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Setting index and positions according to FFTW convention  */
  for( i=0; i<GV.NGRID; i++ ){
    /* REMEMBER kF is (2.0*PI)/L */
    kpos[i]     = (i<GV.NGRID/2) ? GV.KF*i : GV.KF*(i-GV.NGRID);
    indexpos[i] = (i<GV.NGRID/2) ?       i :       (i-GV.NGRID);
  }// for i

  /* These are the minimum and maximum integer vectors 
     possibles according to the grid taken  */
  indexcord[MIN] = -(GV.NGRID/2);
  indexcord[MAX] =  (GV.NGRID/2 - 1);
  
  q1 = (struct densityContrast *) calloc( floor(3 * M_PI * GV.NGRID * GV.NGRID * GV.S_KF), 
  sizeof(struct densityContrast) );
  if(q1 == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: q1 array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ///////////////////////////
  // BISPECTRUM ESTIMATION //
  ///////////////////////////
  
  Nbins = (int) ceil( GV.KN / GV.DELTA_K );
  
  if(rank == 0)
    {
      printf("\n-----------------------------------------------\n");
      printf(" Domain decomposition... \n");
      printf(" %d Total bins\n", Nbins); fflush(stdout);
    }
  
  /* domain decomposition */
  // building bins array
  // declara arreglo que contenga la info de los dk
  // otro inicializado a cero que va a contener el valor de Bk
  bindata = (struct binStruct *) calloc( Nbins, sizeof(struct binStruct) );
  if(bindata == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: bindata array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  //Numero de bines por task = Numero total de bine/size;
  //binspertask = Nbins/size;
  // dar la lista de bines que va a calcular 
  
  taskBin = (int *) calloc( Nbins, sizeof(int) );
  if(taskBin == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: taskBin array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  

  i=0;
  for(l=Nbins-1; l>=0; l--)
    {
      taskBin[l] = i%size;
      
      bindata[l].k1 = (l+0.5)*GV.DELTA_K;
      bindata[l].Nk1 = 0;
      bindata[l].Pk1 = 0.0;         
      bindata[l].Pk1_Error = 0.0;   
      bindata[l].I[RE] = 0.0; 
      bindata[l].I[IM] = 0.0; 
      bindata[l].Ntri = 0;       
      bindata[l].Bk[RE] = 0.0;
      bindata[l].Bk[IM] = 0.0;
      bindata[l].Bk_Error = 0.0; 
      bindata[l].Qk = 0.0;       
      bindata[l].Qk_Error = 0.0;
      bindata[l].Bk_shotnoise = 0.0;
      
      if(rank == 0){
	printf("rank %d has bin %d (k=%f)\n", 
	       taskBin[l], l, bindata[l].k1); fflush(stdout);
      }
      
      i++;
    }
  
  if(rank==0){
    printf("\n-----------------------------------------------\n");
    printf("Bispectrum calculation\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  for(l=0; l<Nbins; l++){
    if(rank != taskBin[l])
      continue;
    
    printf("-----------------------------------------------\nrank:%3d, k1 = %lf\n", 
	   rank, bindata[l].k1);
    fflush(stdout);
    
    //bindata[l].Pk1 = 0.0;
    //bindata[l].Nk1 = 0L;
    for(i=0; i<GV.NGRID; i++){
      for(j=0; j<GV.NGRID; j++){
	for(k=0; k<GV.NGRID; k++){
	  
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	  
	  if( ( bindata[l].k1-GV.DELTA_K*0.5 < kMag ) && 
	      ( kMag < bindata[l].k1+GV.DELTA_K*0.5 ) ){
	    
	    id_cell = INDEX(i,j,k);
	    
	    q1[bindata[l].Nk1].id   = id_cell;
	    q1[bindata[l].Nk1].kMag = kMag;

	    q1[bindata[l].Nk1].triplex[X] = indexpos[i];
	    q1[bindata[l].Nk1].triplex[Y] = indexpos[j];
	    q1[bindata[l].Nk1].triplex[Z] = indexpos[k];
	    
	    bindata[l].Nk1 += 1L;
  
	    bindata[l].Pk1 += COMPLEXMAG(denConK, id_cell);
    
	  }// if bindata[l].k1-GV.DELTA_K*0.5 < kMag < bindata[l].k1+GV.DELTA_K*0.5 ) ){
	  
	}// for k
      }// for j
    }// for i
    
    bindata[l].Pk1 *= (1.0 / bindata[l].Nk1);
    bindata[l].Pk1 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
    bindata[l].Pk1 -= GV.SHOT_NOISE;
    
    /* Estimating bispectrum  */
    //bindata[l].I[RE] = 0.0;
    //bindata[l].I[IM] = 0.0;
    //bindata[l].Ntri  = 0L;
    
    //printf("Elementos en Nk1 %ld\n", bindata[l].Nk1);
    
    for(   rand_i=0;        rand_i<bindata[l].Nk1-1; rand_i++){
      for( rand_j=rand_i+1; rand_j<bindata[l].Nk1;   rand_j++){
	
	/*
	  do{
	  rand_i = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_j = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_k = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  }while( (rand_i==rand_j) || (rand_i==rand_k) || (rand_j==rand_k) );
	*/
	
	/* Looking for closed triangles */
	
	m3[X] = - q1[rand_i].triplex[X] - q1[rand_j].triplex[X];
	m3[Y] = - q1[rand_i].triplex[Y] - q1[rand_j].triplex[Y];
	m3[Z] = - q1[rand_i].triplex[Z] - q1[rand_j].triplex[Z];
	
	if( (indexcord[MIN]<=m3[X] && m3[X]<=indexcord[MAX]) && 
	    (indexcord[MIN]<=m3[Y] && m3[Y]<=indexcord[MAX]) && 
	    (indexcord[MIN]<=m3[Z] && m3[Z]<=indexcord[MAX]) ){
	  
	  i = (m3[X]>=0) ? m3[X] : GV.NGRID+m3[X];
	  j = (m3[Y]>=0) ? m3[Y] : GV.NGRID+m3[Y];
	  k = (m3[Z]>=0) ? m3[Z] : GV.NGRID+m3[Z];
	  
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	    
	  if( ( bindata[l].k1-GV.DELTA_K*0.5 < kMag ) && 
	      ( kMag < bindata[l].k1+GV.DELTA_K*0.5 ) ){
	        
	    id_cell = INDEX(i,j,k);
	    
	    bindata[l].I[RE] += (+ denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][0] * denConK[id_cell][0]
				 - denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][1] * denConK[id_cell][1]
				 - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][0] * denConK[id_cell][1]
				 - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][1] * denConK[id_cell][0] );
	        
	    bindata[l].I[IM] += (+ denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][0] * denConK[id_cell][1]
				 + denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][1] * denConK[id_cell][0]
				 + denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][0] * denConK[id_cell][0]
				 - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][1] * denConK[id_cell][1] );
	        
	    bindata[l].Ntri  += 1L;
	      
	  }// if bindata[l].k1-GV.DELTA_K*0.5 < kMag < bindata[l].k1+GV.DELTA_K*0.5

	}// if

	
	//}// for rand_k
      }// for rand_j
    }// for rand_i
    
    //printf("Numero de triangulos: %ld \n", bindata[l].Ntri);
    
    //if(Ntri==0){ // if Ntri==0, then Bk becomes in a NaN
    //k1 += GV.KF;
    //continue;
    //}
    
    // Integral estimation
    bindata[l].I[RE] *= (1.0 / bindata[l].Ntri);
    bindata[l].I[IM] *= (1.0 / bindata[l].Ntri);

    //printf("Valor de I real: %lf\n", bindata[l].I[RE]);

    // Discrete bispectrum value
    bindata[l].Bk[RE]  = bindata[l].I[RE];
    bindata[l].Bk[RE] *= (GV.SIM_VOL/(1.0*GV.NGRID3)); 
    bindata[l].Bk[RE] *= (GV.SIM_VOL/(1.0*GV.NGRID3));
    bindata[l].Bk[RE] *= (1.0/(1.0*GV.NGRID3));
    
    bindata[l].Bk[IM]  = bindata[l].I[IM];
    bindata[l].Bk[IM] *= (GV.SIM_VOL/(1.0*GV.NGRID3)); 
    bindata[l].Bk[IM] *= (GV.SIM_VOL/(1.0*GV.NGRID3)); 
    bindata[l].Bk[IM] *= (1.0/(1.0*GV.NGRID3));
    
    //printf("Valor del bispectrum con shotnoise: %lf\n", bindata[l].Bk[RE]);
    
    // Stimating shot noise for bispectrum
    bindata[l].Bk_shotnoise = (GV.SHOT_NOISE * 3.0*bindata[l].Pk1) + (GV.SHOT_NOISE*GV.SHOT_NOISE);
    
    /* Substracting shotnoise term  */
    bindata[l].Bk[RE] -= bindata[l].Bk_shotnoise;

    //printf("Valor del bispectrum sin shotnoise: %lf\n", bindata[l].Bk[RE]);
    
    /* Estimating dimensionless bispectrum */
    bindata[l].Qk = bindata[l].Bk[RE] / ( 3.0*(bindata[l].Pk1*bindata[l].Pk1) );

    //Pk1_Error = sqrt( (GV.KF*GV.KF*GV.KF)/(4.0*M_PI*k1*k1*GV.DELTA_K) ) * Pk1;
    //Bk_Error  = sqrt( (GV.KF*GV.KF*GV.KF)/(4.0*M_PI*k1*k1*GV.DELTA_K) ) * Pk1;


    /* Printing bispectrum data  */
    //fprintf(stdout,"%20lf %20e %20e %20e %20e %20ld %20ld %20e %20e\n",
    //bindata[l].k1, bindata[l].Pk1, bindata[l].Bk[RE], bindata[l].Bk[IM], bindata[l].Qk, 
    //bindata[l].Ntri, bindata[l].Nk1, bindata[l].I[RE], bindata[l].I[IM]);
    //fflush(stdout);
    
  }// for l
  //fclose(fout);
  
  
  if(rank==0){
    printf("\n");
    printf("Sending results to root process\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank == 0){

    for(l=1; l<size; l++){
	
      for(j=0; j<Nbins; j++){

	if(l == taskBin[j]){
	  MPI_Recv(&bindata[j], sizeof(struct binStruct), MPI_BYTE, taskBin[j], j, MPI_COMM_WORLD, &status);
	  printf("Receiving bin %d from %d to 0\n",j, rank); fflush(stdout);
	}// if
      }// for j
    }// for l
  }// if rank==0
  else{
    
    for(l=1; l<size; l++){
      
      if(rank == l){
	for(j=0; j<Nbins; j++){
	  if(l == taskBin[j]){
	    printf("Sending bin %d from rank %d to 0\n",j, rank); fflush(stdout);
	    MPI_Send(&bindata[j], (int) sizeof(struct binStruct), MPI_BYTE, 0, j, MPI_COMM_WORLD);
	  }// l==taskBin[j]
	}// for j
      }// rank==l
    }// for l
  }// else
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(rank == 0){
    
    /*
      sprintf(buff,"TodoElResultadoPaNicolas.%d",rank);
      fout = fopen(buff,"w");
      
      for(l=0; l<Nbins; l++)
      fprintf(fout,"%lf %e %e %e %e %ld %ld %e %e\n",
      bindata[l].k1, bindata[l].Pk1, bindata[l].Bk[RE], bindata[l].Bk[IM], bindata[l].Qk, 
      bindata[l].Ntri, bindata[l].Nk1, bindata[l].I[RE], bindata[l].I[IM]);
      
      fclose(fout);
      
      }
    */
  
    /* Saving data in the outfile */
    printf("\n-----------------------------------------------\n");
    printf("Saving data in %s\n", GV.OUTPUT);
  
    
    fout = fopen(GV.OUTPUT, "w");
    if(fout == NULL){
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Outfile could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
    
    /* Writing header */
    fprintf(fout,"# NGRID          = %d\n",  GV.NGRID);
    fprintf(fout,"# GADGET VERSION = %d\n",  GV.GADGET_VERSION);
    fprintf(fout,"# L              = %lf\n", GV.L);
    fprintf(fout,"# SIM VOL        = %lf\n", GV.SIM_VOL);
    fprintf(fout,"# NP TOT         = %ld\n", GV.NP_TOT);
    fprintf(fout,"# TOTAL MASS     = %lf\n", GV.TOTAL_MASS);
    fprintf(fout,"# RHO MEAN       = %lf\n", GV.RHO_MEAN);
    fprintf(fout,"# VOL_CELL       = %lf\n", GV.VOL_CELL);
    fprintf(fout,"# H              = %lf\n", GV.H);
    fprintf(fout,"# DELTA k        = %lf\n", GV.DELTA_K);
    fprintf(fout,"# kF             = %lf\n", GV.KF);
    fprintf(fout,"# kN             = %lf\n", GV.KN);
    fprintf(fout,"# Shot Noise     = %lf\n", GV.SHOT_NOISE);
    fprintf(fout,"# SCHEME         = %s\n",  GV.SCHEME);
    fprintf(fout,"# OMEGA_M0       = %lf\n", GV.OMEGA_M0);
    fprintf(fout,"# OMEGA_L0       = %lf\n", GV.OMEGA_L0);
    fprintf(fout,"# ZRS            = %lf\n", GV.ZRS);
    fprintf(fout,"# HUBBLEPARAM    = %lf\n", GV.HUBBLEPARAM);
    fprintf(fout,"\n");
    
    fprintf(fout,"#%19s %20s %20s %20s %20s %20s %20s %20s %20s\n",
	    "k", "P(k)", "Re[B(k)]", "Im[B(k)]", "Q(k)", 
	    "Ntri", "Nk", "Re[I]", "Im[I]");
    
    /* Printing bispectrum data  */
    for(l=0; l<Nbins; l++){
      fprintf(fout,"%20lf %20e %20e %20e %20e %20ld %20ld %20e %20e\n",
	      bindata[l].k1, bindata[l].Pk1, bindata[l].Bk[RE], bindata[l].Bk[IM], bindata[l].Qk, 
	      bindata[l].Ntri, bindata[l].Nk1, bindata[l].I[RE], bindata[l].I[IM]);
    }// for l
    fclose(fout);
    
  }// if rank==0
    
  
  
  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
  //fclose(fout);
  free(kpos);
  free(indexpos);
  fftw_free(denConK);
  free(q1);
  free(bindata);
  free(taskBin);
  
  
  
  /////////////////////
  //* FINISHING MPI *//
  /////////////////////
  MPI_Finalize();
  
  
  
  return 0;

}

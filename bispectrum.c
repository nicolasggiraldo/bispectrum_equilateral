#include     <stdio.h>
#include    <stdlib.h>
#include      <math.h>
#include <libconfig.h>
#include     <fftw3.h>
#include    <string.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <unistd.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

/*
  FFTW_MEASURE: find the optimal plan by actually computing several FFTs 
  FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
  FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
  FFTW_IN_PLACE: a plan assumes that the in and out are same
*/

int main(int argc, char *argv[]){
  int i, j, k;              // Counter in the X, Y, and Z axis.
  
  unsigned long int rand_i; /* Index for the subsample of density contrast */
  unsigned long int rand_j; /* in the range between                        */
  //unsigned long int rand_k; /* k1-DELTA_K/2 to k1+DELTA_K/2                */
  long int id_cell;         // Counter for the id of the cell
  double k1;                /* Side of the triangle to evaluate the 
			       Bispectrum equilateral configuration */
  struct densityContrast *q1=NULL; /* Array of structure with the subset of 
				      density contrast with range between 
				      k1-DELTA_K/2 to k1+DELTA_K/2 */
  long int Nk1;             // Number of elements in array q1
  long int Ntri;            // Number of counted triangles for bispectrum estimation
  double Pk1;               // Power spectrum at k1
  double I[2];              /* Sum of elements (denConk1*denConk1*denConk1) 
			       which form a triangle. Because denConk1 is a 
			       complex value the sum is also complex */
  double Bk[2];             /* Bispectrum value to measure from the N-body 
			       snapshot in general the bispectrum is a complex
			       quantity but as the avarage of density constrast
			       is taken the imaginary part goes to zero (see 
			       Gil-Marin et al. 2012, J. Cosm. Astrop. Phys) */
  double Qk;                /* Reduced bispectrum, this value will be evaluated
			       with the real part of the bispectrum */
  double Bk_shotnoise;      /* Bispectrum shotnoise. Bispectrum shot noise emerge 
			       as an effect of discretization of positions */
  double kMag;              /* Magnitude of the wavenumber vector */
  FILE *fout = NULL;        // File handler to output
  double *kpos;             /* Array of positions values according to FFTW 
			       k-position convention */ 
  int *indexpos;            /* Array of index values according to FFTW 
			       k-position convention */
  int indexcord[2];
  int m3[3];
  long int posTri;          // Number of posible triangles

  /*
  gsl_rng *r=NULL;
  long seed;
  const gsl_rng_type *T;

  seed = time(NULL) * getpid();
  T = gsl_rng_gfsr4;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  //*/
  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("%s: You must specify the name of the parameter file\n",argv[0]);
    printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1]) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad path to the parameter file.\n" );
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
  switch( readBinaryFile() ){
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

  /* Do forward FFT */
  fftw_execute(forwardPlan);
  printf("\n-----------------------------------------------\n");
  printf("Fourier Transform succes.\n");
  fftw_free(denConX);

  /* Position array for storing in the densityContrast  */
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
  indexcord[MIN] = -(GV.NGRID/2);
  indexcord[MAX] =  (GV.NGRID/2 - 1);
  
  q1 = (struct densityContrast *) calloc( GV.NGRID3, sizeof(struct densityContrast) );
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
  printf("\n-----------------------------------------------\n");
  printf("Bispectrum calculation\n");

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
	  "k", "P(k)", "Re[B(k)]", "Im[B(k)]", "Q(k)", "Ntri", "Nk", "Re[I]", "Im[I]");

  /* Lineal binning */
  k1 = GV.DELTA_K;
  
  //while(k3 <= GV.KN){
  while( k1 <= GV.KN ){
    
    printf("k1 = %lf\n", k1);
    fflush(stdout);
    
    Pk1 = 0.0;
    Nk1 = 0L;
    for(i=0; i<GV.NGRID; i++){
      for(j=0; j<GV.NGRID; j++){
	for(k=0; k<GV.NGRID; k++){
	  
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	  id_cell = INDEX(i,j,k);
	  
	  if( ( k1-GV.DELTA_K*0.5 < kMag ) && ( kMag < k1+GV.DELTA_K*0.5 ) ){
	    
	    q1[Nk1].id   = id_cell;
	    q1[Nk1].kMag = kMag;

	    q1[Nk1].triplex[X] = indexpos[i];
	    q1[Nk1].triplex[Y] = indexpos[j];
	    q1[Nk1].triplex[Z] = indexpos[k];
	    
	    Nk1 += 1L;
	  
	    Pk1 += COMPLEXMAG(denConK, id_cell);
	    
	  }// if
	  
	}// for k
      }// for j
    }// for i
    
    Pk1 *= (1.0/Nk1);
    Pk1 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
    Pk1 -= GV.SHOT_NOISE;
    
    /* Estimating bispectrum  */
    I[RE] = 0.0;
    I[IM] = 0.0;
    Ntri  = 0L;
    
    posTri = ( Nk1 * (Nk1-1) * (Nk1-2) ) / 6;
    printf("Elementos en Nk1 %ld\nPosibles triangulos=%ld\n", Nk1, posTri);
    
    //for(id_cell=0L; id_cell<posTri; id_cell+=1L){
    for(   rand_i=0;        rand_i<Nk1-1; rand_i++){
      for( rand_j=rand_i+1; rand_j<Nk1;   rand_j++){
	
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
	  
	  if( ( k1-GV.DELTA_K*0.5 < kMag ) && ( kMag < k1+GV.DELTA_K*0.5 ) ){
	    
	    id_cell = INDEX(i,j,k);
	    
	    //if( ( (q1[rand_i].triplex[X] + q1[rand_j].triplex[X] + q1[rand_k].triplex[X]) == 0 ) &&
	    //    ( (q1[rand_i].triplex[Y] + q1[rand_j].triplex[Y] + q1[rand_k].triplex[Y]) == 0 ) &&
	    //    ( (q1[rand_i].triplex[Z] + q1[rand_j].triplex[Z] + q1[rand_k].triplex[Z]) == 0 ) ){
	
	    I[RE] += ( + denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][0] * denConK[id_cell][0]
		       - denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][1] * denConK[id_cell][1]
		       - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][0] * denConK[id_cell][1]
		       - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][1] * denConK[id_cell][0] );
	    
	    I[IM] += ( + denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][0] * denConK[id_cell][1]
		       + denConK[q1[rand_i].id][0] * denConK[q1[rand_j].id][1] * denConK[id_cell][0]
		       + denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][0] * denConK[id_cell][0]
		       - denConK[q1[rand_i].id][1] * denConK[q1[rand_j].id][1] * denConK[id_cell][1] );
	    
	    Ntri    += 1L;
	  
	  }// if
	}// if

	
	//}// for rand_k
      }// for rand_j
    }// for rand_i
    //}// for id_cell
    
    printf("Valor de I real despues del ciclo: %lf", I[RE]);
    printf("Numero de triangulos: %ld \n", Ntri);
    
    if(Ntri==0){ // if Ntri==0, then Bk becomes in a NaN
      k1 += GV.KF;
      continue;
    }
    
    // Integral estimation by Monte Carlo
    I[RE] *= (1.0/Ntri);
    I[IM] *= (1.0/Ntri);
    //I *= (GV.DELTA_K*GV.DELTA_K*GV.DELTA_K);

    printf("Valor de I real despues de multiplicar: %lf", I[RE]);

    // Discrete bispectrum value
    Bk[RE] = I[RE] * (GV.SIM_VOL/(1.0*GV.NGRID3)) * (GV.SIM_VOL/(1.0*GV.NGRID3)) * (1.0/(1.0*GV.NGRID3));
    Bk[IM] = I[IM] * (GV.SIM_VOL/(1.0*GV.NGRID3)) * (GV.SIM_VOL/(1.0*GV.NGRID3)) * (1.0/(1.0*GV.NGRID3));
    
    printf("Valor del bispectrum con shotnoise: %lf\n", Bk[RE]);
    
    // Stimating shot noise for bispectrum
    Bk_shotnoise = (GV.SHOT_NOISE * 3.0*Pk1) + (GV.SHOT_NOISE*GV.SHOT_NOISE);
    
    /* Substracting shotnoise term  */
    Bk[RE] -= Bk_shotnoise;

    printf("Valor del bispectrum sin shotnoise: %lf\n", Bk[RE]);

    /* Estimating dimensionless bispectrum */
    Qk = Bk[RE] / ( 3.0*(Pk1*Pk1) );

    /* Printing bispectrum data  */
    fprintf(fout,"%20lf %20e %20e %20e %20e %20ld %20ld %20e %20e\n",
	    k1, Pk1, Bk[RE], Bk[IM], Qk, Ntri, Nk1, I[RE], I[IM]);
    fflush(fout);
    
    k1 += GV.KF;
    
  }// while
  printf("\n");
  
  
  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
  fclose(fout);
  free(kpos);
  free(indexpos);
  fftw_free(denConK);
  free(q1);
  fftw_destroy_plan(forwardPlan);
  
  return 0;
}

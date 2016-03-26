/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded are: 
 * FILE_NAME:      File name path of the GADGET binary file.
 * OUTPUT:         Path of the output file.
 * DELTA_K:        Width of space sampled for the calculation of 
 *                 the power spectrum. The value is given in terms 
 *                 of the fundamental frequency kF. The value
 *                 should be bigger or equal than 1.
 * 
 * The parameter file is read with the help of the library libconfig, 
 * for more information of the library use go to the manual:
 * http://www.hyperrealm.com/libconfig/libconfig_manual.html
 *
 *
 *  param_file_name: String with the name of the parameter file.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[]){

  config_t cfg;            /* Returns all parameters in this structure */
  const char *str1, *str2; /* Going to be used to read strings variables */
  

  /*Initialization */
  config_init(&cfg);
  
  /* Read the file. If there is an error, report it and exit. */
  if(!config_read_file(&cfg, param_file_name)){
    printf("%s:%d - %s\n",
	   config_error_file(&cfg),
	   config_error_line(&cfg),
	   config_error_text(&cfg));
    config_destroy(&cfg);
    // Value -1 means there is an error loading the param file
    return -1;
  }

  /* Get the value of the width sample in terms of the 
     fundamental frequency kF. */
  if( config_lookup_float(&cfg, "S_KF", &(GV.S_KF) ) ){
    // Checking if NGRID value is valid, that is, NGRID > 0
    if(GV.S_KF >= 1.0){
      printf("Binning width in terms of the fundamental frequency kF: %lf\n", GV.S_KF);
    }
    else{
      printf("Invalid 'S_KF' setting in configuration file.\n");
      return -2;
    }
  }
  else{
    printf("No 'S_KF' setting in configuration file.\n");
    return -2;
  }
  
  /* Get the configuration file name. */
  if(config_lookup_string(&cfg, "FILE_NAME", &str1)){
    GV.FILE_NAME = strdup(str1);
    printf("Reading from File: %s\n", GV.FILE_NAME);
  }
  else{
    printf("No 'FILE_NAME' setting in configuration file.\n");
    return -2;
  }
  
  /* Get the configuration output. */
  if(config_lookup_string(&cfg, "OUTPUT", &str2)){
    GV.OUTPUT = strdup(str2);
    printf("Output File: %s\n", GV.OUTPUT);
  }
  else{
    printf("No 'OUTPUT' setting in configuration file.\n");
    return -2;
  }
  
  
  config_destroy(&cfg);
  
  
  return 0;
}



/*                                                                                                                      
 * Function:  readBinaryFile                                                                                      
 * --------------------                                                                                                 
 * Reads a binary file with the information of the cell and store 
 * the information in the data structure variable *part* it also 
 * returns the total number of particles.
 *                                                                                                                      
 *  There are no arguments in the routiene.                                                                             
 *                                                                                                                      
 *  returns: Integer value.                                                                                             
 *            0 --> There is no error.                                                                                  
 *           -1 --> There is an error loading file
 *           -2 --> Structure cell could not be allocated.    
 */
int readBinaryFile(){
  FILE *fdata;
  long int id_cell;
  double density_Contrast;
  int n[3]; // Number of grids in each axis for FFT estimation
  size_t err;

  printf("\n-----------------------------------------------\n");
  printf("Reading file:   %s\n", GV.FILE_NAME);
  
  fdata = fopen(GV.FILE_NAME,"rb");
  if(fdata == NULL){
    printf("File %s cannot be open\n", GV.FILE_NAME);
    return -1;
  }


  /* Getting cosmological parameters of the simulation */
  err = fread(&GV.OMEGA_M0,    sizeof(double), 1, fdata);
  err = fread(&GV.OMEGA_L0,    sizeof(double), 1, fdata);
  err = fread(&GV.ZRS,         sizeof(double), 1, fdata);
  err = fread(&GV.HUBBLEPARAM, sizeof(double), 1, fdata);

  /* Getting simulation parameters */
  err = fread(&GV.NGRID,          sizeof(int),      1, fdata);
  err = fread(&GV.GADGET_VERSION, sizeof(int),      1, fdata);
  err = fread(&GV.L,              sizeof(double),   1, fdata);
  err = fread(&GV.NP_TOT,         sizeof(long int), 1, fdata);
  err = fread(&GV.TOTAL_MASS,     sizeof(double),   1, fdata);
  err = fread(&GV.RHO_MEAN,       sizeof(double),   1, fdata);
  err = fread(&GV.VOL_CELL,       sizeof(double),   1, fdata);
  err = fread(&GV.H,              sizeof(double),   1, fdata);
  err = fread(&(GV.SCHEME[0]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[1]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[2]),    sizeof(char),     1, fdata);
  GV.SCHEME[3] = '\0';

  GV.NGRID3     = (1L*GV.NGRID) * (1L*GV.NGRID) * (1L*GV.NGRID);
  GV.SIM_VOL    = GV.L * GV.L * GV.L;
  GV.KF         = (2.0*M_PI) / GV.L;
  GV.DELTA_K    = GV.S_KF * GV.KF;
  GV.SHOT_NOISE = GV.VOL_CELL / GV.NP_TOT;
  GV.KN         = M_PI / GV.H;

  printf("\n-----------------------------------------------\n");
  printf("The original snapshot has a total of %ld particles\n", GV.NP_TOT);
  printf("----------------------------------------\n");
  printf(" * Redshift...     %16.8lf\n", GV.ZRS);
  printf(" * Omega0...       %16.8lf\n", GV.OMEGA_M0);
  printf(" * OmageLa...      %16.8lf\n", GV.OMEGA_L0);
  printf(" * Hubbleparam...  %16.8lf\n", GV.HUBBLEPARAM);
  printf("----------------------------------------\n");
  printf(" * Boxsize...      %16.8lf\n", GV.L);
  printf(" * Ngrid...        %16d\n",    GV.NGRID);
  printf(" * SimMass...      %16.8e\n", GV.TOTAL_MASS);
  printf(" * Scheme...       %16s\n",    GV.SCHEME);
  printf("----------------------------------------\n");
  printf(" * kF...           %16.8lf\n", GV.KF);
  printf(" * kN...           %16.8lf\n", GV.KN);
  printf(" * DELTA_k...      %16.8lf\n", GV.DELTA_K);
  printf(" * PSshotNoise...  %16.8e\n", GV.SHOT_NOISE);

  /* Memory allocation for the input and output FFTW arrays,
     for the space, initialize input array to (1.,0.) */
  denConX = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  denConK = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConX == NULL || denConK == NULL){
    printf("FFTW arrays could not be allocated\n");
    return -2;
  }//if

  /************************/
  /* FFTW plan allocation */
  /************************/
  
  /* Number of grids in each axis */
  n[X] = n[Y] = n[Z] = GV.NGRID;
  
  /* FFTW plan for a 3D Fourier transform */
  forwardPlan = fftw_plan_dft(3, n, denConX, denConK,
			      FFTW_FORWARD, FFTW_ESTIMATE);

  /****************************/
  /* Getting density contrast */
  /****************************/
  for(id_cell=0L; id_cell<GV.NGRID3; id_cell++){
    err = fread(&density_Contrast, sizeof(double), 1, fdata);
    denConX[id_cell][0] = density_Contrast;
    denConX[id_cell][1] = 0.0;
  }//for id_cell

  if(err){};

  return 0;
}

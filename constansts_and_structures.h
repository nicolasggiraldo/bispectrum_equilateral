/////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTTURES FOR THE PROGRAM //
/////////////////////////////////////////////////////



/////////////////////////////
// PREPROCESSOR DIRECTIVES //
/////////////////////////////

/* Index preprocessor for the C-Order */
#define INDEX(i,j,k)     (1L*k) + GV.NGRID * ( (1L*j) + GV.NGRID * (1L*i) )
#define COMPLEXMAG(A,i)  ( (A[i][0] * A[i][0]) + (A[i][1] * A[i][1]) )
#define VECTORMAG(x,y,z) sqrt( (x * x) + (y * y) + (z * z) )

#define RE 0 // Preprocessor directive to indicate the real part of a complex value
#define IM 1 // Preprocessor directive to indicate the imaginary part of a complex value

#define MIN 0 // Preprocessor directive to indicate the minimum value of a set
#define MAX 1 // Preprocessor directive to indicate the maximum value of a set

#define X  0
#define Y  1
#define Z  2



////////////////////////////////////////////
// GLOBAL VARIABLES FOR THE FFTW ROUTINES //
////////////////////////////////////////////

// FFTW variables
fftw_complex *denConX = NULL; // Density contrast in X-space
fftw_complex *denConK = NULL; // Density contrast in K-space
fftw_plan forwardPlan;



///////////////////////////////
// STRUCTURES OF THE PROGRAM //
///////////////////////////////

/* Density contrast */
struct densityContrast{
  long int id;
  double kMag;
  int triplex[3];
};

/* Global variales */
struct globalVariables{
  int      NGRID;          // Number of cell in each axis.
  long int NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  int      GADGET_VERSION; // GADGET version of the snapshot
  double   L;              // Lenght of the simulation in Mpc
  double   SIM_VOL;        // Volume of the simulation
  long int NP_TOT;         // Total number of particles in the simulation
  double   TOTAL_MASS;     // Total mass of all particles in the simulation
  double   RHO_MEAN;       // Mean density of ALL the simulation
  double   VOL_CELL;       // Volume of each cell
  double   H;              // Size of the cell
  double   S_KF;           // Constant to define the value of DELTA_K = S_KF * KF
  double   DELTA_K;        // Delta k for the binning
  double   KF;             // Fundamental frequency kF
  double   SHOT_NOISE;     // shotNoise of the Power Spectrum
  double   KN;             // Nyquist frequency
  char     *FILE_NAME;     // Path of the GADGET binary
  char     *OUTPUT;        // Name of the outputfile
  char     SCHEME[4];      // Scheme used for grid assignation

  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  double OMEGA_M0;         //Omega matter at present time
  double OMEGA_L0;         //Omega Lambda at present time
  double ZRS;              //Redshift of the simulation
  double HUBBLEPARAM;      //Hubble parameter of the simulation
  
}GV;

/* data for output */
struct binStruct{
  double k1;           /* Side of the triangle to evaluate the                                                     
			  Bispectrum equilateral configuration */
  long int Nk1;        // Number of elements in array q1 
  double Pk1;          // Power spectrum at k1 
  double Pk1_Error;    // Power spectrum error at k1 
  double I[2];         /* Sum of elements denConk1*denConk2*denConk3 which 
			  form a triangle. Because denConk is a complex 
			  quantity the sum is also complex */
  long int Ntri;       // Number of counted triangles for bispectrum estimation  
  double Bk[2];        /* Bispectrum value to measure from the simulation 
			  snapshot in general the bispectrum is a complex                                   
			  quantity but as the avarage of density constrast                           
			  is taken the imaginary part goes to zero (see                       
			  Gil-Marin et al. 2012, J. Cosm. Astrop. Phys) */
  double Bk_Error;     // Bispectrum error at k1, k2, k3
  double Qk;           /* Reduced bispectrum, this value will be evaluated
			  with the real part of the bispectrum only */
  double Qk_Error;     // Reduced bispectrum error
  double Bk_shotnoise; // Bispectrum shotnoise
};

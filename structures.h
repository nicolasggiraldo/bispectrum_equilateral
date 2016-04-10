////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTURES FOR THE PROGRAM //
////////////////////////////////////////////////////

/* Global variales */
struct globalVariables{
  int NGRID;           // Number of cell in one axis.
  double L;            // Lenght of the simulation in Mpc.
  double SimVol;       // Volume of the simulation
  int NpTot;           // Total number of particles
  int NGRID3;          // Total number of cells
  double rhoMean;      // Mean density of ALL the simulation
  double mass;         // Mass of each particle
  double volCell;      // Volume of each cell
  double dx;           // Also known as H : Size of the cell
  char FILENAME[1000]; // Path of the gadget file
  char OUTPUT[1000];   // Name of the outputfile
  double deltaK;       // Delta k for the binning
  double kF;
  double shotNoise;    // shotNoise of the Power Spectrum
  double kN;           // Nyquist frequency

  /* COSMOLOGICAL PARAMETERS*/
  double OmegaM0;     //Omega matter at present time
  double OmegaL0;     //Omega Lambda at present time
  double zRS;         //Redshift of the simulation
  double HubbleParam; //Hubble parameter of the simulation

}GV;

/*
struct Cell{
  int       Np_cell; // Number of particles in the Cell.
  long int idCell; // Array with the ID of the particles inside the cell.
  double       posx;
  double       posy;
  double       posz;
  double       velx; // Mean velocity of the particles in the X-Axis 
  double       vely; // Mean velocity of the particles in the Y-Axis 
  double       velz; // Mean velocity of the particles in the Z-Axis 
  double     denCon; // Density contrast of the cell
  double        rho; // Density in the cell
}*cells;
*/

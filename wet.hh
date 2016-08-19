double sign(const double value);                                                                // looks ok

void LGConfig();

void nextToDo();
double ECRate = 0.1;
int restartPoint = 1;

void marker();
double *posX, *posY, *posZ;
int markers;
int modes = 0;

void initialiseSurface();
int readInput();
void LBAlgorithm(int runStep);
void saveFiles();

void savePartialDensity();                                                                       // written 20/1/2005 for Parallel
void loadPartialDensity();                                                                       // written 12/3/2005 for Parallel
void generateDrops();                                                                            // written 16/6/2005
void evaporate();
void initialise();                                                                               // checked
void setBodyForce();                                                                             // looks ok

void equilibrium();                                                                              // checked
void collision();                                                                                // checked
void propagation();                                                                              // looks ok
void applyBoundaryConditions();                                                                  // rewritten for 2-D flow, need to change for 3-D flow
void synchronisefn();										 // written on 20/1/2005 for 3-D simulation
void synchroniseff();										 // written on 20/1/2005 for 3-D simulation

void computeCoordinate();                                                                        // looks ok
void computeMomenta();                                                                           // looks ok
void computeFreeEnergy();                                                                        // checked
double getTotalDensity();                                                                        // looks ok
double getDropletDensity();                                                                      // written 17/12/2004 only for serial

void computeCentreOfMass();                                                                      // checked 17/12/2004
void writeDensityFile();                                                                         // looks ok : ddt.m
void writeYPlanDensityFile(const int iy);                                                        // looks ok : dyjtt.m
void writeDensityFileOnSubstrate();                                                              // looks ok : st.m
void writeYPlanXVelocityFile(const int iy);                                                      // looks ok : uxYjtt.m and uzYjtt.m

#ifdef PARALLEL                                                                                  // added 31/12/2005
void exchangeBoundaries();
void sendBoundaries();
void receiveBoundaries();
void exchangeDensitiesAtBoundaries();
void sendDensitiesAtBoundaries();
void receiveDensitiesAtBoundaries();
void exchangeVelocitiesAtBoundaries();
void exchangePressureAtBoundaries();
#endif

double beta = 2.0, pc = 0.04, nc = 3.5, Tc = 4.0 / 7.0;                                   // Parameters of the chosen free energy density, defined in initialise()
																						  //const double beta = 0.5, pc = 1.0/2.0, nc = 5.0, Tc = 4.0/7.0;                                   // Parameters of the chosen free energy density, defined in initialise()
const double w0 = 8.0/27.0, w1 = 2.0 / 27.0, w2 = 1.0 / 54.0, w3 = 1/216.0, wa1 = w1, wa2 = w2;  // corresponds to w_\sigma, used in equilibrium(), setBodyForce().
double lw, lw1, lw2, lw3;                                                                        // log of the w_i, used in collision()
double dx = 1.0;                                                                                 // To change this must change derivatives
double dt = 1.0;
double c = dx / dt, c2 = c*c;

int LX, LY, LZ, N;                                                                                  // The number of lattice points, defined in wet.dat
char *mask, *subMask;                                                                             // Lattice ID for topolocially patterned substrate, defined in initialise().
int nbEqStep, PDStep, nbIter, infoStep;                                                             // Number of iterations and steps at which info is collected, defined in wet.dat
double nGas, nFluid;                                                                              // Densities defined in wet.dat. Actually: 2.51 - 4.49. Some numbers quoted earlier 2.61 - 4.55
double tau, tauLiquid, tauGas, kappa;                                                                                // tau = relaxation time, kappa = constant related to grad in density, defined in wet.dat
double nu;                                                                                       // Kinematic Viscosity. defined in initialise()
double T;                                                                                        // Temperature, defined in initialise()
double phi11, phi12, phi180;                                                                              // Surface free energy parameters corresponding E'm contact angles, calculated in initialise()
double teta1, teta2;                                                                              // E'm contact angles, defined in wet.dat 
int LPattern1, LPattern2;                                                                        // Pattern length for striped and patch surface
int xwidth, ywidth, zwidth;

String loadedFile;

int dropletR, dropletCenterX, dropletCenterY, dropletCenterZ;                                       // Initial parameters of the droplet, defined in wet.dat
int dropletR2, dropletCenterX2, dropletCenterY2, dropletCenterZ2;                                       // Initial parameters of the droplet, defined in wet.dat
double initUX, initUY, initUZ;                                                                     // Initial velocity of the droplet, defined in wet.dat
double energy;                                                                                   // Total energy of the system, used in computeFreeEnergy()
double freeE;                                                                                    // Bulk Free Energy, used in computeFreeEnergy()
double surfaceE;                                                                                 // Surface Energy
double interfaceE;                                                                               // Interface between Gas and Liquid
double nue;                                                                                      // (n-nc)/nc, used in equilibrium() and computeFreeEnergy()
double tau_w, gam;                                                                                // (Tc-T)/Tc, defined in initialise()
double G[3], GG[15];                                                                              // G the force in 3-direction. GG force as it applies to the partial distribution function

																								  // These are probably obsolete, not sure though 
long *dd1, *dd2, *dd3, *dd4, *dd5, *dd6, *dd7, *dd8, *dd9, *dda, *ddb, *ddc, *ddd, *dde;         // The neighbours of each lattice point, -1 boundary, defined in initialise()
long d1, d2, d3, d4, d5, d6, d7, d8, d9, da, db, dc, dd, de;                                     // Temporary variables of dd
double *ff0, *ff1, *ff2, *ff3, *ff4, *ff5, *ff6, *ff7, *ff8, *ff9, *ffa, *ffb, *ffc, *ffd, *ffe; // After streaming
double *fn0, *fn1, *fn2, *fn3, *fn4, *fn5, *fn6, *fn7, *fn8, *fn9, *fna, *fnb, *fnc, *fnd, *fne; // After collision, before streaming
double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, fa, fb, fc, fd, fe;                               // Temporary variables of ff
double fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8, fe9, fea, feb, fec, fed, fee;                // Equilibrium partial density

int vl[15][3];                                                                                   // Velocity direction, defined in initialise()
long *dd1, *dd2, *dd3, *dd4, *dd5, *dd6, *dd7, *dd8, *dd9, *dda, *ddb, *ddc, *ddd, *dde;         // The neighbours of each lattice point, -1 boundary, defined in initialise()
long *dd0, *dda1, *dda2, *dda3, *dda4, *dda5, *dda6, *ddb1, *ddb2, *ddb3, *ddb4, *ddb5, *ddb6;         // ASA: The neighbours of each lattice point, -1 boundary, defined in initialise()
long *ddb7, *ddb8, *ddb9, *ddba, *ddbb, *ddbc, *ddc1, *ddc2, *ddc3, *ddc4, *ddc5, *ddc6, *ddc7, *ddc8;
long d1, d2, d3, d4, d5, d6, d7, d8, d9, da, db, dc, dd, de;                                     // Temporary variables of dd
double *ff0, *ffa1, *ffa2, *ffa3, *ffa4, *ffa5, *ffa6, *ffb1, *ffb2, *ffb3, *ffb4, *ffb5, *ffb6; // ASA: After streaming
double *ffb7, *ffb8, *ffb9, *ffba, *ffbb, *ffbc, *ffc1, *ffc2, *ffc3, *ffc4, *ffc5, *ffc6, *ffc7, *ffc8;
double *fn0, *fna1, *fna2, *fna3, *fna4, *fna5, *fna6, *fnb1, *fnb2, *fnb3, *fnb4, *fnb5, *fnb6; // ASA: After collision, before streaming
double *fnb7, *fnb8, *fnb9, *fnba, *fnbb, *fnbc, *fnc1, *fnc2, *fnc3, *fnc4, *fnc5, *fnc6, *fnc7, *fnc8;
double f0, fa1, fa2, fa3, fa4, fa5, fa6, fb1, fb2, fb3, fb4, fb5, fb6;                           // ASA: Temporary variables of ff
double fb7, fb8, fb9, fba, fbb, fbc, fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8;
double fe0, fea1, fea2, fea3, fea4, fea5, fea6, feb1, feb2, feb3, feb4, feb5, feb6;              // ASA: Equilibrium partial density
double feb7, feb8, feb9, feba, febb, febc, fec1, fec2, fec3, fec4, fec5, fec6, fec7, fec8;
double dfe0, dfea1, dfea2, dfea3, dfea4, dfea5, dfea6, dfeb1, dfeb2, dfeb3, dfeb4, dfeb5, dfeb6; // ASA: Difference between euilibrium and regular partial density
double dfeb7, dfeb8, dfeb9, dfeba, dfebb, dfebc, dfec1, dfec2, dfec3, dfec4, dfec5, dfec6, dfec7, dfec8;
double tfe0, tfea1, tfea2, tfea3, tfea4, tfea5, tfea6, tfeb1, tfeb2, tfeb3, tfeb4, tfeb5, tfeb6; // ASA: Temporary variables for fe
double tfeb7, tfeb8, tfeb9, tfeba, tfebb, tfebc, tfec1, tfec2, tfec3, tfec4, tfec5, tfec6, tfec7, tfec8;
double *n, *uxs, *uys, *uzs;                                                                     // Density and velocity of lattice points
double ux, uy, uz, ux2, uy2, uz2, u2, uxy, uxz, uyz;                                                            // Temporary velocity variables 
double del2n;                                                                                    // Laplacian, used in equilibrium()
double nn;                                                                                       // Temporary density variable
double dn_dx, dn_dy, dn_dz;                                                                      // Density gradient
double comX = 0, comY = 0, comZ = 0;                                                             // Centre of mass positions in all 3 directions
double comVX, comVY, comVZ;                                                                      // Centre of mass velocities in all 3 directions
int crossing;                                                                                    // How many times the droplet has passed the periodic boundary condition
																								 // ASA's variables 
long d11, d12, d13, d14, d15, d16, d21, d22, d23, d24, d25, d26, d31, d32, d33, d34, d35, d36,   // Neighbours of neighbours of each lattice point
d41, d42, d43, d44, d45, d46, d51, d52, d53, d54, d55, d56, d61, d62, d63, d64, d65, d66;
double dxx, dxy, dxz, dyy, dyz, dzz, dthird;                                                             // Second order partial derivatives of density n
double dpr_dn;                                                                                   // derivative of pressure w.r.t density
double alpha;                                                                                    // Alpha from the paper: term to modify collision
double cs2, cs2_div, cs4_div, cs6_div;                                                           // Various multiples of cs
double apr, bpr, rg;                                                                             // Constants in EoS
double T2, T3, T4, T5, T6;                                                                       // Powers of temperature
double ka;                                                                                       // constant k in pressure tensor in paper- not sure what value to give it
double al0, al1, alm;                                                                            // Variables for alpha in Eq. 4 of paper, and alm is 1 - alpha/tau
int h;                                                                                           // dummy index in Newton-Rahpson method
double entropy;                                                                                  // Entropy of f + alpha*(feq - f) minus entropy of f
double dentropy;                                                                                 // Derivative of the above

double a;                                                                                        // A_\sigma / w_\sigma
double a1, a2, a3;                                                                               // corresponds to G_{1\gamma\gamma} v_{i\gamma} v_{i\gamma}
double a4, a5, a6, a7;                                                                           // corresponds to G_{2\gamma\delta} v_{i\gamma} v_{i\delta}
double cxy, cyz, czx;                                                                              // corresponds to \kappa(d_\gamma n)(d_\delta n)+\nu(u_\gamma d_\delta n + u_\delta d_\gamma n)

int xk, yk, zk;                                                                                    // 3-D index version of index k
int myPE, nbPE;                                                                                   // Number of processors and processor ID
long k1, k2, k11, k22;                                                                              // Start and End index of lattice sites
long i, j, h;                                                                                    // Running dummy index in x, y, and z direction
long k;                                                                                          // Index of current lattice point
long t;                                                                                          // Running dummy index for number of iterations
long d, u, r, l, q, w;                                                                           // dummy index that works out the neighbours of k {(l,r),(d,u),(q,w)}

double *nGlobal;                                                                                 // Used in writeDensityFile() 

char *buff;
int i2;

/*
double*uxGlobal,*uzGlobal,*pxxGlobal,*pxzGlobal;
double del_sq_n;
double b, pp, del_sq_p;
double dp_dx, dp_dy, dp_dz;                                                                      // Pressure gradient
double *pxx,*pxz,*pzz;

const double A = 9.0/49, B = 2.0/21;
char *buff;
int **postsHeight;
double bodyForce;
double s_lg,s_sl,s_sg;
double patchDensity;    // Density of the posts or patches
int width1,width2;      // Width and Length of the posts
*/

#ifdef PARALLEL                                                                                  // Added 31/12/2005
MPI_Status status;
MPI_Datatype leftFieldsType, rightFieldsType;
int leftNeighbor, rightNeighbor;
#endif


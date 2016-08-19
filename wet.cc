#ifdef PARALLEL
#include <mpi.h>
#undef	SEEK_SET
#undef	SEEK_CUR
#undef	SEEK_END
#endif

#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <sys/time.h>

#include "RT_Timer.hh"
#include "String.hh"

#include "wet.hh"


void marker()
{
  computeMomenta();
  int tempk = 0;
  for (int a = 0; a < markers; a++) {
    tempk = int(posX[a] - int(posX[a]/LX)*LX)*LZ*LY + int(posY[a] - int(posY[a]/LY)*LY)*LZ + int(posZ[a] - int(posZ[a]/LZ)*LZ);
    posX[a] += uxs[tempk];
    posY[a] = 0; //+= uys[tempk];
    posZ[a] += uzs[tempk];

    if(t % infoStep == 0){
    std::cout<<tempk << " " << posX[a] << " " <<uxs[tempk] <<std::endl;
    }
  }
        
  if(t % infoStep == 0){
    String fileName("markers.dat");
    std::ofstream fout(fileName.get(), std::ios::app);
    fout.precision(4);
    for (int a = 0; a < markers; a++)
      {
	fout<<posX[a]<<" "<<posY[a]<<" "<<posZ[a]<<" ";
      }
    fout<<std::endl;
    fout.close();
  }


  /*int tempk = 0;
  int tempx = 0; int tempy = 0; int tempz = 0;
  double ddx = 0.0; double ddy = 0.0; double ddz = 0.0;
  double dux_dx, duy_dx, duz_dx,dux_dy, duy_dy, duz_dy,dux_dz, duy_dz, duz_dz; 
  int a = 0; 
  
  for (a = 0; a < markers; a++) {

    tempx = int(posX[a]+0.5); tempy = int(posY[a]+0.5); tempz = int(posZ[a]+0.5);
    tempk = int(tempx)*LZ*LY + int(tempy)*LZ + int(tempz);
    ddx = posX[a] - tempx; ddy = posY[a] - tempy; ddz = posZ[a] - tempz;
    d1 = dd1[tempk]; d2 = dd2[tempk]; d3 = dd3[tempk]; d4 = dd4[tempk]; d5 = dd5[tempk]; d6 = dd6[tempk];

    dux_dx = (uxs[d1] - uxs[d2])/2;
    duy_dx = (uys[d1] - uys[d2])/2;
    duz_dx = (uzs[d1] - uzs[d2])/2;
    dux_dy = (uxs[d3] - uxs[d4])/2;
    duy_dy = (uys[d3] - uys[d4])/2;
    duz_dy = (uzs[d3] - uzs[d4])/2;
    dux_dz = (uxs[d5] - uxs[d6])/2;
    duy_dz = (uys[d5] - uys[d6])/2;
    duz_dz = (uzs[d5] - uzs[d6])/2;

    posX[a] += uxs[tempk] + ddx*dux_dx + ddy*dux_dy + ddz*dux_dz;
    posY[a] += uys[tempk] + ddx*duy_dx + ddy*duy_dy + ddz*duy_dz;
    posZ[a] += uzs[tempk] + ddx*duz_dx + ddy*duz_dy + ddz*duz_dz;
    if(posX[a] < 0) posX[a] += LX;
    if(posX[a] >= LX) posX[a] -= LX;
    if(posY[a] < 0) posY[a] += LY;
    if(posY[a] >= LY) posY[a] -= LY;
    if(posZ[a] < 0) posZ[a] += LZ;
    if(posZ[a] > LZ) posZ[a] -= LZ;
  }

  if(t % infoStep == 0){
    String fileName("markers.dat");
    std::ofstream fout(fileName.get(), std::ios::app);
    fout.precision(4);
    for (a = 0; a < markers; a++)
    {
	fout<<posX[a]<<" "<<posY[a]<<" "<<posZ[a]<<" ";
    }
    fout<<std::endl;
    fout.close();
    }*/

}


// =============================================================
/* This function informs the code what to do next             */
// =============================================================

void nextToDo()
{

  std::ifstream inputFile("ToDo.dat");
  if (!inputFile) {
    std::cout << "Can't open the input file" << std::endl;
  }
  String endOfLine;
  String ContinueOrNot;
  int iterationSteps;

  inputFile >> restartPoint >> endOfLine;
  inputFile >> ECRate >> endOfLine;
  inputFile >> iterationSteps >> endOfLine;
  inputFile >> ContinueOrNot >> endOfLine;

  if (ContinueOrNot == "No") {
    t = restartPoint;
    nbIter = iterationSteps;
    loadPartialDensity();
  }

}

// =============================================================
/* This function initialise the Liquid gas configuration      */
// =============================================================

#ifdef SPHERICAL
void LGConfig() 
{
    
  /*double a = (xk-dropletCenterX)*(xk-dropletCenterX) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-dropletCenterZ)*(zk-dropletCenterZ);
  double b = (xk-dropletCenterX2)*(xk-dropletCenterX2) + (yk-dropletCenterY2)*(yk-dropletCenterY2) + (zk-dropletCenterZ2)*(zk-dropletCenterZ2);

  if ( (sqrt(a) <= dropletR || sqrt(b) <= dropletR2) ) { 
    n[k] = nFluid; 
    uxs[k] = initUX; uys[k] = initUY; uzs[k] = initUZ;
  }
  else {
    n[k]=nGas;
    uxs[k] = 0.0; uys[k] = 0.0; uzs[k] = 0.0;   
    }*/

  double a = (xk-dropletCenterX)*(xk-dropletCenterX) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-dropletCenterZ)*(zk-dropletCenterZ);
  if (sqrt(a) <= dropletR) {
    n[k] = nFluid; 
    uxs[k] = initUX; uys[k] = initUY; uzs[k] = initUZ;
  }
  else {
    n[k]=nGas;
    uxs[k] = 0.0; uys[k] = 0.0; uzs[k] = 0.0;	  
  }
}
#endif

#ifdef SQUARE
void LGConfig() 
{
    int random1 = 0, random2 = 0;
    if ( xk%(LPattern1+LPattern2) > LPattern1 ) random1 = int(xk/(LPattern1+LPattern2)) * 0;
    if ( xk%(LPattern1+LPattern2) > LPattern1 ) random2 = abs(int(xk/(LPattern1+LPattern2)-1) * 0);

    if ( abs(xk-dropletCenterX) <= xwidth && (yk-dropletCenterY) <= (ywidth-random1) && (yk-dropletCenterY) >= (random2-ywidth) && abs(zk-dropletCenterZ) <= zwidth ) {
      n[k] = nFluid;
      uxs[k] = initUX; uys[k] = initUY; uzs[k] = initUZ;
    }
    else {
      n[k]=nGas;
      uxs[k] = 0.0; uys[k] = 0.0; uzs[k] = 0.0;	  
   }
}
#endif

// =============================================================
/* This function defines chemically striped surface           */
// =============================================================
#ifdef STRIPED

void initialiseSurface()
{

    if (zk == 0) 
    {
    	mask[k]=26;

	if (yk > -1) {
	  if ( xk%(LPattern1+LPattern2) < LPattern1 )
	    subMask[k] = 0;
	  else
	    subMask[k] = 1;
	}
	else subMask[k] = 2;

    } 
    if (zk == LZ-1) 
    {
    	mask[k]=27;
	
	if ( xk%(LPattern1+LPattern2) < LPattern1 )
	  subMask[k] = 0;
	else
	  subMask[k] = 1;
    }   
}

#endif

#ifdef PATCH

void initialiseSurface()
{

      if (zk == 0) 
    {
    	mask[k]=26;

        if ( (xk%(LPattern1+LPattern2) < LPattern1) && (yk%(LPattern1+LPattern2) < LPattern1) )
           subMask[k] = 0;
        else
          subMask[k] = 1;

    } 
    if (zk == LZ-1) 
    {
    	mask[k]=27;
	
        if ( (xk%(LPattern1+LPattern2) < LPattern1) && (yk%(LPattern1+LPattern2) < LPattern1) )
           subMask[k] = 0;
        else
          subMask[k] = 1;

	  }   
}

#endif

// =============================================================
/* This function defines the drop separator geometry          */
// =============================================================
#ifdef SEPARATOR

void initialiseSurface()
{
    double VStripe = 20;
    double HStripe1 = 20; double HStripe2 = 30; double HStripe3 = 40;
    double ISSV = 200; double ISSH = 80;

    if (zk == 0) 
    {
    	mask[k]=26;
	if ( xk < ISSV )
	  subMask[k] = 0;
	else if ( xk >= ISSV && xk < (ISSV+VStripe) )
	  subMask[k] = 1;
	else if ( xk >= (ISSV+VStripe) && xk < (2*ISSV+VStripe) )
	  subMask[k] = 0;
	else if ( xk >= (2*ISSV+VStripe) && xk < (2*ISSV+2*VStripe) ) 
	  subMask[k] = 1;
	else
	subMask[k] = 0;

	if ( yk >= ISSH && yk < (ISSH+HStripe1) )
	  subMask[k] = 1;
	if ( yk >= (2*ISSH+HStripe1) && yk < (2*ISSH+HStripe1+HStripe2) ) 
	  subMask[k] = 1;
	if ( yk >= (3*ISSH+HStripe1+HStripe2) && yk < (3*ISSH+HStripe1+HStripe2+HStripe3) ) 
	  subMask[k] = 1;
    }

    if (zk == LZ-1) 
    {
    	mask[k]=27;
	if ( xk < ISSV )
	  subMask[k] = 0;
	else if ( xk >= ISSV && xk < (ISSV+VStripe) )
	  subMask[k] = 1;
	else if ( xk >= (ISSV+VStripe) && xk < (2*ISSV+VStripe) )
	  subMask[k] = 0;
	else if ( xk >= (2*ISSV+VStripe) && xk < (2*ISSV+2*VStripe) ) 
	  subMask[k] = 1;
	else
	subMask[k] = 0;

	if ( yk >= ISSH && yk < (ISSH+HStripe1) )
	  subMask[k] = 1;
	if ( yk >= (2*ISSH+HStripe1) && yk < (2*ISSH+HStripe1+HStripe2) ) 
	  subMask[k] = 1;
	if ( yk >= (3*ISSH+HStripe1+HStripe2) && yk < (3*ISSH+HStripe1+HStripe2+HStripe3) ) 
	  subMask[k] = 1;
    } 
}

#endif

// =============================================================
/* This function defines the drop generator geometry          */
// =============================================================
#ifdef GENERATOR

void initialiseSurface()
{
    if (zk == 0) 
    {
    	mask[k]=26;
	if ( xk < 60 )
	  subMask[k] = 0;
	else if ( xk >= 60 && xk < 130)
	  subMask[k] = 1;
	else if ( xk >= 130 && xk < 160)
	  subMask[k] = 0;
	else
	subMask[k] = 1;
    } 
    if (zk == LZ-1) 
    {
    	mask[k]=27;
	if ( xk < 60 )
	  subMask[k] = 0;
	else if ( xk >= 60 && xk < 130)
	  subMask[k] = 1;
	else if ( xk >= 130 && xk < 160)
	  subMask[k] = 0;
	else
	subMask[k] = 1;
    } 
}

#endif

// =============================================================
/* This main program investigates 3-D liquid drop flow on     */
/* patterned substrate                                        */
// =============================================================

int main(int argc, char **argv)
{

// ##########################################################
/* Define non-physical parameters of the simulation        */
// ##########################################################

#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
#else
  nbPE=1;
  myPE=0;
#endif

// ##########################################################
/* Initialisation and print initial information            */
// ##########################################################

  readInput();
  std::cout<<"Read Input"<<std::endl;
  initialise();
  std::cout<<"Initialise"<<std::endl;

  if (myPE == 0) {
  	std::ofstream file("energy.m", std::ios::app);
	file.close();
	std::ofstream fout("density.m",std::ios::app);
  	fout.precision(4);
	fout.close();
  	std::cout << "Start: Number of processors = " << nbPE << " - Initial density = " << getTotalDensity() << std::endl;
  }
  else {
  	getTotalDensity();
  }
   
// #########################################################
/* Spreading Steps                                        */
// #########################################################

  if (loadedFile == "No") {
	
    t = 0;
    writeYPlanDensityFile(dropletCenterY);
    writeDensityFileOnSubstrate();
    writeDensityFile();
    computeCentreOfMass();

    double Gtemp1 = G[0], Gtemp2 = G[1], Gtemp3 = G[2];
    G[0] = 0.0; G[1] = 0.0; G[2] = 0.0;
    setBodyForce();

    modes = 0;
    LBAlgorithm(nbEqStep);
    
/* set the G[0] back                                     */
    G[0] = Gtemp1; G[1] = Gtemp2; G[2] = Gtemp3;
    setBodyForce();
    savePartialDensity();

  }
  else {
    nextToDo();
    std::cout << "ReadToDo" << std::endl;
    loadPartialDensity();
    std::cout << "Loading ... " << loadedFile.get() << std::endl;
  }

// #########################################################
/* Flowing Steps                                          */
// #########################################################

  markers = 20;
  posX = new double [markers];
  posY = new double [markers];
  posZ = new double [markers];
  for (int a = 0; a < markers; a++) {
    posX[a] = dropletCenterX;
    posY[a] = dropletCenterY;
    posZ[a] = 5 + 2*a;
  }

  modes = 1;
  LBAlgorithm(nbIter);

#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
  
}


// =============================================================
/* This function reads all the input parameters at wet.dat    */
// =============================================================

int readInput()
{

  std::ifstream inputFile("wet.dat");
  if (!inputFile) {
    std::cout << "Can't open the input file" << std::endl;
    return 1;
  }
  String endOfLine;
  inputFile >> LX >> endOfLine;
  inputFile >> LY >> endOfLine;
  inputFile >> LZ >> endOfLine;
  inputFile >> nbEqStep >> endOfLine;
  inputFile >> PDStep >> endOfLine;
  inputFile >> nbIter >> endOfLine;
  inputFile >> infoStep >> endOfLine;
  inputFile >> tauLiquid >> endOfLine;
  inputFile >> tauGas >> endOfLine;
  inputFile >> kappa >> endOfLine;
  inputFile >> T >> endOfLine;
  inputFile >> beta >> endOfLine;
  inputFile >> pc >> endOfLine;
  nFluid = nc * (1 + sqrt(beta*0.3));
  nGas   = nc * (1 - sqrt(beta*0.3));
  inputFile >> G[0] >> G[1] >> G[2] >> endOfLine;
  inputFile >> loadedFile >> endOfLine;
  inputFile.close();

#ifdef STRIPED
   
    std::ifstream inputFile2("striped.dat");
    if (!inputFile2) {
        std::cout << "Can't open the input file" << std::endl;
    }
    inputFile2 >> LPattern1 >> endOfLine;
    inputFile2 >> LPattern2 >> endOfLine;
    inputFile2 >> teta1 >> endOfLine;
    inputFile2 >> teta2 >> endOfLine;
    teta1=teta1/180.0*M_PI; // Conversion from deg to rad
    teta2=teta2/180.0*M_PI; // Conversion from deg to rad

#endif

#ifdef PATCH
   
    std::ifstream inputFile2("striped.dat");
    if (!inputFile2) {
        std::cout << "Can't open the input file" << std::endl;
    }
    inputFile2 >> LPattern1 >> endOfLine;
    inputFile2 >> LPattern2 >> endOfLine;
    inputFile2 >> teta1 >> endOfLine;
    inputFile2 >> teta2 >> endOfLine;
    teta1=teta1/180.0*M_PI; // Conversion from deg to rad
    teta2=teta2/180.0*M_PI; // Conversion from deg to rad

#endif

#ifdef SPHERICAL

    std::ifstream inputFile3("spherical.dat");
    if (!inputFile3) {
        std::cout << "Can't open the input file" << std::endl;
    }
    inputFile3 >> dropletR >> endOfLine; // use -1 for no droplet
    inputFile3 >> dropletCenterX >> endOfLine;
    inputFile3 >> dropletCenterY >> endOfLine; // always use 1 in 3D simulation using symmetry
    inputFile3 >> dropletCenterZ >> endOfLine;
    inputFile3 >> initUX >> initUY >> initUZ >> endOfLine;

    inputFile3 >> dropletR2 >> endOfLine; // use -1 for no droplet
    inputFile3 >> dropletCenterX2 >> endOfLine;
    inputFile3 >> dropletCenterY2 >> endOfLine; // always use 1 in 3D simulation using symmetry
    inputFile3 >> dropletCenterZ2 >> endOfLine;

#endif

#ifdef SQUARE

    std::ifstream inputFile3("square.dat");
    if (!inputFile3) {
        std::cout << "Can't open the input file" << std::endl;
    }
    inputFile3 >> xwidth >> ywidth >> zwidth >> endOfLine; // use -1 for no droplet
    inputFile3 >> dropletCenterX >> endOfLine;
    inputFile3 >> dropletCenterY >> endOfLine; // always use 1 in 3D simulation using symmetry
    inputFile3 >> dropletCenterZ >> endOfLine;
    inputFile3 >> initUX >> initUY >> initUZ >> endOfLine;
    dropletR = int(pow(3*xwidth*ywidth*zwidth/4/M_PI,1/3));

#endif

    return 0;
} 


// =============================================================
/* This function runs the spreading/flowing steps             */
// =============================================================

void LBAlgorithm(int runStep)
{

    for(t = restartPoint; t <= runStep; t++){ 

      computeMomenta();
#ifdef PARALLEL
      exchangeDensitiesAtBoundaries();
#endif
      collision();    
#ifdef SYMMETRY
      synchronisefn();
#endif
#ifdef PARALLEL
      exchangeBoundaries();
#endif
      propagation();
      applyBoundaryConditions();
#ifdef SYMMETRY
      synchroniseff();
#endif
      saveFiles();
      if (modes == 1) marker();

#ifdef VOLUME 
      if(t % PDStep == 0) {
	    generateDrops();
      }
#endif      

      if(t % infoStep == 0) {
	nextToDo();
      }

    }

}

// ============================================================
/* This function saves the required files                    */
// ============================================================

void saveFiles()
{
      if(t % infoStep == 0){
	computeFreeEnergy();
	if (myPE == 0) {
	  std::ofstream file("energy.m",std::ios::app);
	  file.precision(12);
	  file << energy << " " << freeE<< " " << interfaceE << " " << surfaceE << std::endl;
	  file.close();
	  std::ofstream fout("density.m",std::ios::app);
	  fout.precision(4);
	  fout << t << " " << getDropletDensity() << std::endl;
	  fout.close();
	}
	else {
	  getDropletDensity();
	}
	writeYPlanDensityFile(dropletCenterY);
	writeYPlanXVelocityFile(dropletCenterY);
	if (LY > 1) {
	  writeDensityFileOnSubstrate();
	  writeDensityFile();
	}
	computeCentreOfMass();
      }
      
      if(t % PDStep == 0){
      	savePartialDensity();
      }

}

// ============================================================
/* This function returns the sign of some input value        */
// ============================================================

double sign(const double value)
{
  if (value < 0.0) return -1.0;
  return 1.0;
}

// ===============================================================
/* This function saves the partial distribution function        */
// ===============================================================


void savePartialDensity()
{
  
  String fileName("PD1P");
  fileName.concat((int)myPE);
  fileName.concat("t");
  fileName.concat((int) t);
  fileName.concat(".m");
  std::ofstream fout(fileName.get());
  fout.precision(8);

  for (k = 0; k < N; k++)
  {
    fout<<ff0[k]<<" "<<ff1[k]<<" "<<ff2[k]<<" "<<ff3[k]<<" "<<ff4[k]<<" "<<ff5[k]<<" "<<ff6[k]<<" "<<ff7[k]<<" "<<ff8[k]<<" "<<ff9[k]<<" "<<ffa[k]<<" "<<ffb[k]<<" "<<ffc[k]<<" "<<ffd[k]<<" "<<ffe[k]<<std::endl;
  }

  fout.close();

}


// ===============================================================
/* This function loads the partial distribution function        */
// ===============================================================

void loadPartialDensity()
{

  String endOfLine;
  String fileName("PD1P");
  fileName.concat((int)myPE);
  fileName.concat("t");
  fileName.concat((int)restartPoint);
  fileName.concat(".m");
  std::ifstream fin(fileName.get());
  
  for (k = 0; k < N; k++)
    {
      fin>>ff0[k]>>ff1[k]>>ff2[k]>>ff3[k]>>ff4[k]>>ff5[k]>>ff6[k]>>ff7[k]>>ff8[k]>>ff9[k]>>ffa[k]>>ffb[k]>>ffc[k]>>ffd[k]>>ffe[k]>>endOfLine;
    }

  fin.close();

}

// ===============================================================
/* This function generates small drops        */
// ===============================================================

void generateDrops()
{

  for (k=0;k<N;k++) {

    computeCoordinate();
    
    /*if (zk >= LZ-15) {
      n[k] += 0.2;
      nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
      ff0[k] = fe0; 
      ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
      ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
      ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
      ffd[k] = fed; ffe[k] = fee; 
      }*/

    /*double a = (xk-dropletCenterX)*(xk-dropletCenterX) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-dropletCenterZ)*(zk-dropletCenterZ);
    if (sqrt(a) <= (dropletR-5) && zk > postHeight) {
    	n[k] += 0.2;
    	nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
    	ff0[k] = fe0; 
    	ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
    	ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
    	ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
    	ffd[k] = fed; ffe[k] = fee; 
	}*/

    if (n[k] >= nc) {
    	n[k] += ECRate;
    	nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
    	ff0[k] = fe0; 
    	ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
    	ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
    	ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
    	ffd[k] = fed; ffe[k] = fee; 
    }

  }   

}

// ===============================================================
/* This function evaporates the drop        */
// ===============================================================

void evaporate()
{  
  for (k=0;k<N;k++) {

    computeCoordinate();
    
    /*if (zk >= LZ-15) {
      n[k] -= 0.2;
      nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
      ff0[k] = fe0; 
      ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
      ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
      ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
      ffd[k] = fed; ffe[k] = fee; 
      }*/

    /*double a = (xk-dropletCenterX)*(xk-dropletCenterX) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-dropletCenterZ)*(zk-dropletCenterZ);
    if (sqrt(a) <= (dropletR-5) && zk > postHeight) {
    	n[k] -= 0.2;
    	nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
    	ff0[k] = fe0; 
    	ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
    	ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
    	ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
    	ffd[k] = fed; ffe[k] = fee; 
	}*/

    if (n[k] >= nc) {
    	n[k] -= ECRate;
    	nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; equilibrium();
    	ff0[k] = fe0; 
    	ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
    	ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
    	ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
    	ffd[k] = fed; ffe[k] = fee; 
    }

  }   

}

// ===============================================================
/* This function initialises 2-D liquid drop flow on chemically */
/* patterned substrate under constant force in the x-direction  */
// ===============================================================

void initialise()
{  

// ##########################################################
/* Define the limit of the loop variables                  */
// ##########################################################

  int xoffset=0;
#ifdef PARALLEL
  N=(LX/nbPE+2)*LY*LZ;
  k1=LZ*LY;
  k2=N-LZ*LY;
  k11=k1+LZ*LY;
  k22=k2-LZ*LY;
  xoffset=1;    
  nGlobal=new double[LX*LY*LZ];
  i2=LX/nbPE+2;
#else
  N=LX*LY*LZ;
  k1=0;
  k2=N;
  k11=0;
  k22=0;
  i2=LX;
#endif

// ##########################################################
/* Define the physical constants                           */
// ##########################################################

  tau_w=(Tc-T)/Tc;
  gam= 1.0-(T/Tc);

// ##########################################################
/* Allocate memory for variables                           */
// ##########################################################

  ff0 = new double[N]; ff1 = new double[N]; ff2= new double[N];
  ff3 = new double[N]; ff4 = new double[N]; ff5= new double[N];
  ff6 = new double[N]; ff7 = new double[N]; ff8= new double[N];
  ff9 = new double[N]; ffa = new double[N]; ffb= new double[N];
  ffc = new double[N]; ffd = new double[N]; ffe= new double[N];

  fn0 = new double[N]; fn1 = new double[N]; fn2= new double[N];
  fn3 = new double[N]; fn4 = new double[N]; fn5= new double[N];
  fn6 = new double[N]; fn7 = new double[N]; fn8= new double[N];
  fn9 = new double[N]; fna = new double[N]; fnb= new double[N];
  fnc = new double[N]; fnd = new double[N]; fne= new double[N];

  n = new double[N];

  uxs = new double[N];
  uys = new double[N];
  uzs = new double[N];

  dd1 = new long[N]; dd2 = new long[N]; dd3 = new long[N]; 
  dd4 = new long[N]; dd5 = new long[N]; dd6 = new long[N]; 
  dd7 = new long[N]; dd8 = new long[N]; dd9 = new long[N]; 
  dda = new long[N]; ddb = new long[N]; ddc = new long[N]; 
  ddd = new long[N]; dde = new long[N]; 

  mask = new char[N];
  subMask = new char[N];

// ###########################################################
/* Define the mask and submask for each lattice point.      */
/* Initialise density and velocity. Initial 3D droplet      */
/* located on (dropletCenterX,dropletCenterY,dropletCenterZ)*/ 
/* with a radius equal to dropletR.                         */
/* Find the system boundary and denote as -1.               */
// ###########################################################

  for (k=0;k<N;k++) {

    mask[k] = 0;
    subMask[k] = 1;

    computeCoordinate();

    LGConfig();
    initialiseSurface();
// ###########################################################
/* Define the neighbours.                                   */
// ###########################################################

    int xkLattice = int(k/(float) (LZ*LY)); // xkLattice = xk for serial program
    
#ifdef PARALLEL
    l = xkLattice-1; r = xkLattice+1;
#else
    if (xkLattice == 0) l = LX - 1;
    else l = xkLattice - 1; 
    if (xkLattice == LX -1) r = 0;
    else r = xkLattice + 1; 
#endif
    
    if (yk == 0) d = LY - 1;
    else d = yk - 1; 
    if (yk == LY -1) u = 0;
    else u = yk + 1; 

    q = zk-1; w = zk+1;  // will apply appropriate boundary condition later on          

    /*    if (zk == 0) q = LZ - 1;
    else q = zk-1; 
    if (zk == LZ-1) w = 0;
    else w = zk+1;  */
    
    dd1[k] = zk + yk*LZ + r*LZ*LY; 	  
    dd2[k] = zk + yk*LZ + l*LZ*LY; 
    dd3[k] = zk + u*LZ + xkLattice*LZ*LY; 
    dd4[k] = zk + d*LZ + xkLattice*LZ*LY; 
    dd5[k] = w + yk*LZ + xkLattice*LZ*LY; 
    dd6[k] = q + yk*LZ + xkLattice*LZ*LY; 
    dd7[k] = w + u*LZ + r*LZ*LY; 
    dd8[k] = q + d*LZ + l*LZ*LY; 
    dd9[k] = w + u*LZ + l*LZ*LY; 
    dda[k] = q + d*LZ + r*LZ*LY; 
    ddb[k] = w + d*LZ + l*LZ*LY; 
    ddc[k] = q + u*LZ + r*LZ*LY; 
    ddd[k] = w + d*LZ + r*LZ*LY; 
    dde[k] = q + u*LZ + l*LZ*LY; 


    if (zk == 0)
      {
	dd6[k] = -1;
	dd8[k] = -1;
	dda[k] = -1;
	ddc[k] = -1;
	dde[k] = -1;
      }
    
    if (zk == LZ - 1)
      {
	dd5[k] = -1;
	dd7[k] = -1;
	dd9[k] = -1;
	ddb[k] = -1;
	ddd[k] = -1;
      } 

#ifdef PARALLEL
    if (xkLattice == 0)
    {
	dd2[k] = -1;
	dd8[k] = -1;
	dd9[k] = -1;
	ddb[k] = -1;
	dde[k] = -1;	
    }    
    if (xkLattice == i2 - 1)
    {
	dd1[k] = -1;
	dd7[k] = -1;
	dda[k] = -1;
	ddc[k] = -1;
	ddd[k] = -1;    
    }    
#endif

  }

// ###########################################################
/* Initialise the partial distribution function             */
// ###########################################################

  dn_dx = 0; dn_dy = 0; dn_dz = 0; del2n = 0;
  
  for(k = k1; k < k2; k++){
  
    if (mask[k] == 28) n[k]=1.0;
    nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; 
    equilibrium();
    ff0[k] = fe0; 
    ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; 
    ff5[k] = fe5; ff6[k] = fe6; ff7[k] = fe7; ff8[k] = fe8; 
    ff9[k] = fe9; ffa[k] = fea; ffb[k] = feb; ffc[k] = fec;
    ffd[k] = fed; ffe[k] = fee; 
  
  }

// ###########################################################
/* Velocity vector definition                               */
// ###########################################################

  vl[0][0] = 0;     vl[0][1]= 0;     vl[0][2]= 0; 

  vl[1][0] = 1;     vl[1][1]= 0;     vl[1][2]= 0; 
  vl[2][0] =-1;     vl[2][1]= 0;     vl[2][2]= 0; 
  vl[3][0] = 0;     vl[3][1]= 1;     vl[3][2]= 0; 
  vl[4][0] = 0;     vl[4][1]=-1;     vl[4][2]= 0; 
  vl[5][0] = 0;     vl[5][1]= 0;     vl[5][2]= 1; 
  vl[6][0] = 0;     vl[6][1]= 0;     vl[6][2]=-1;
 
  vl[7][0] = 1;     vl[7][1]= 1;     vl[7][2]= 1; 
  vl[8][0] =-1;     vl[8][1]=-1;     vl[8][2]=-1; 
  vl[9][0] =-1;     vl[9][1]= 1;     vl[9][2]= 1; 
  vl[10][0]= 1;    vl[10][1]=-1;    vl[10][2]=-1;
  vl[11][0]=-1;    vl[11][1]=-1;    vl[11][2]= 1;
  vl[12][0]= 1;    vl[12][1]= 1;    vl[12][2]=-1;
  vl[13][0]= 1;    vl[13][1]=-1;    vl[13][2]= 1;
  vl[14][0]=-1;    vl[14][1]= 1;    vl[14][2]=-1;

// ###########################################################
/* Set the body force             */
// ###########################################################

  setBodyForce();

// ###########################################################
/* Calculate the surface free energy constants              */
/* Here phi11 = - phi1 / kappa (in literature)              */
/* And  phi12 = - phi2 / kappa (in literature)              */
// ###########################################################

  double alpha=acos(sin(teta1)*sin(teta1));
  phi11=2.0*beta*tau_w*sqrt(2.0*pc/kappa)*sign(teta1-M_PI/2.0)*
       sqrt(cos(alpha/3.0)*(1.0-cos(alpha/3.0)));

  alpha=acos(sin(teta2)*sin(teta2));
  phi12=2.0*beta*tau_w*sqrt(2.0*pc/kappa)*sign(teta2-M_PI/2.0)*
        sqrt(cos(alpha/3.0)*(1.0-cos(alpha/3.0)));

  alpha=acos(sin(M_PI)*sin(M_PI));
  phi180=2.0*beta*tau_w*sqrt(2.0*pc/kappa)*sign(M_PI-M_PI/2.0)*
        sqrt(cos(alpha/3.0)*(1.0-cos(alpha/3.0)));

#ifdef PARALLEL

  // Creates a communication buffer
  const long int bufSize=10000000;//LY*LZ*5*8+5000;
  buff=new char[bufSize];
  MPI_Buffer_attach(buff,bufSize);
  // One should deattach it at the end...

  leftNeighbor=myPE-1;
  if (leftNeighbor == -1) leftNeighbor=nbPE-1;

  rightNeighbor=myPE+1;
  if (rightNeighbor == nbPE) rightNeighbor=0;

#endif
  
}

// ============================================================
/* This function set the body force of the system            */
// ============================================================

void setBodyForce()
{

  double fx = G[0]/c;
  double fy = G[1]/c;  
  double fz = G[2]/c;
   
  double fsq = fx*ux+fy*uy+fz*uz;
  double Mxy = fx*uy + fy*ux;
  double Myz = fy*uz + fz*uy;
  double Mzx = fz*ux + fx*uz;

  GG[1] = w1*(3*fx*ux - fsq) + w1*fx;
  GG[2] = w1*(3*fx*ux - fsq) - w1*fx;
  GG[3] = w1*(3*fy*uy - fsq) + w1*fy;
  GG[4] = w1*(3*fy*uy - fsq) - w1*fy;
  GG[5] = w1*(3*fz*uz - fsq) + w1*fz;
  GG[6] = w1*(3*fz*uz - fsq) - w1*fz;
  GG[7] = w2*(2*fsq + 3*( Mxy + Myz + Mzx)) + w2*( fx + fy + fz);
  GG[8] = w2*(2*fsq + 3*( Mxy + Myz + Mzx)) + w2*(-fx - fy - fz);
  GG[9] = w2*(2*fsq + 3*(-Mxy + Myz - Mzx)) + w2*(-fx + fy + fz);
  GG[10] = w2*(2*fsq + 3*(-Mxy + Myz - Mzx)) + w2*( fx - fy - fz);
  GG[11] = w2*(2*fsq + 3*( Mxy - Myz - Mzx)) + w2*(-fx - fy + fz);
  GG[12] = w2*(2*fsq + 3*( Mxy - Myz - Mzx)) + w2*( fx + fy - fz);
  GG[13] = w2*(2*fsq + 3*(-Mxy - Myz + Mzx)) + w2*( fx - fy + fz);
  GG[14] = w2*(2*fsq + 3*(-Mxy - Myz + Mzx)) + w2*(-fx + fy - fz);
  GG[0] = -GG[1]-GG[2]-GG[3]-GG[4]-GG[5]-GG[6]-GG[7]-GG[8]-GG[9]-GG[10]-GG[11]-GG[12]-GG[13]-GG[14];

}

// ================================================================
/* This function calculates the equilibrium partial distribution */
/* for a given value of k                                        */
/*---------------------------------------------------------------*/
/* a = A_\sigma / w_\sigma                                       */
/* a1,a2,a3 = G_{1\gamma\gamma} v_{i\gamma} v_{i\gamma}          */
/* cxy,cyz,cxz = \kappa(d_\gamma n)(d_\delta n) +                */
/*               \nu(u_\gamma d_\delta n + u_\delta d_\gamma n)  */
/* a4,a5,a6,a7 = G_{2\gamma\delta} v_{i\gamma} v_{i\delta}       */
// ================================================================

void equilibrium()
{

  ux2 = ux*ux; uy2 = uy*uy; uz2 = uz*uz; 
  u2=ux2+uy2+uz2;
  double nue= nn/nc - 1.0;
  
  tau = tauGas + (tauLiquid - tauGas)*(1+tanh(3*(nn-nc)))/2;
  nu = dx*c2*(2*tau - 1)/6;
  
  a = ((pc*(nue+1.0)*(nue+1.0)*(3.0*nue*nue-2.0*nue+1.0-2.0*beta*gam))
       - kappa*((dn_dx*dn_dx + dn_dy*dn_dy + dn_dz*dn_dz)/2)
       + nu*(ux*dn_dx + uy*dn_dy + uz*dn_dz))/c2 - kappa*nn*del2n/c2;

  a1 = (kappa*dn_dx*dn_dx+2*nu*ux*dn_dx)/(2*c2);
  a2 = (kappa*dn_dy*dn_dy+2*nu*uy*dn_dy)/(2*c2);
  a3 = (kappa*dn_dz*dn_dz+2*nu*uz*dn_dz)/(2*c2);
  cxy = kappa*dn_dx*dn_dy + nu*(ux*dn_dy + uy*dn_dx);
  cyz = kappa*dn_dy*dn_dz + nu*(uy*dn_dz + uz*dn_dy);
  czx = kappa*dn_dz*dn_dx + nu*(uz*dn_dx + ux*dn_dz);
  a4 = ( cxy + cyz + czx)/(8*c2);
  a5 = (-cxy + cyz - czx)/(8*c2);
  a6 = ( cxy - cyz - czx)/(8*c2);
  a7 = (-cxy - cyz + czx)/(8*c2);

  fe1 = w1*(nn*( ux + 1.5*ux2 - 0.5*u2)) + wa1*a + a1;
  fe2 = w1*(nn*(-ux + 1.5*ux2 - 0.5*u2)) + wa1*a + a1;
  fe3 = w1*(nn*( uy + 1.5*uy2 - 0.5*u2)) + wa1*a + a2;
  fe4 = w1*(nn*(-uy + 1.5*uy2 - 0.5*u2)) + wa1*a + a2;
  fe5 = w1*(nn*( uz + 1.5*uz2 - 0.5*u2)) + wa1*a + a3;
  fe6 = w1*(nn*(-uz + 1.5*uz2 - 0.5*u2)) + wa1*a + a3;
  fe7 = w2*(nn*( ux + uy + uz + u2 + 3*( ux*uy + uy*uz + uz*ux))) + wa2*a + a4;
  fe8 = w2*(nn*(-ux - uy - uz + u2 + 3*( ux*uy + uy*uz + uz*ux))) + wa2*a + a4;
  fe9 = w2*(nn*(-ux + uy + uz + u2 + 3*(-ux*uy + uy*uz - uz*ux))) + wa2*a + a5;
  fea = w2*(nn*( ux - uy - uz + u2 + 3*(-ux*uy + uy*uz - uz*ux))) + wa2*a + a5;
  feb = w2*(nn*(-ux - uy + uz + u2 + 3*( ux*uy - uy*uz - uz*ux))) + wa2*a + a6;
  fec = w2*(nn*( ux + uy - uz + u2 + 3*( ux*uy - uy*uz - uz*ux))) + wa2*a + a6;
  fed = w2*(nn*( ux - uy + uz + u2 + 3*(-ux*uy - uy*uz + uz*ux))) + wa2*a + a7;
  fee = w2*(nn*(-ux + uy - uz + u2 + 3*(-ux*uy - uy*uz + uz*ux))) + wa2*a + a7;
  fe0 = nn - fe1 - fe2 - fe3 - fe4 - fe5 - fe6 - fe7 - fe8 - fe9 - fea - feb 
      - fec - fed - fee;
      
}

// ================================================================
/* This function performs the collision step of LB simulation    */
/* calculation of dn_dx and del2n have been modified from Alex's */
// ================================================================

void collision()
{

  for( k = k1 ; k < k2 ; k++){
    
    if (mask[k] != 28) {

	    computeCoordinate();
	
	    d1 = dd1[k]; d2 = dd2[k]; d3 = dd3[k]; d4 = dd4[k]; d5 = dd5[k]; d6 = dd6[k];
	    d7 = dd7[k]; d8 = dd8[k]; d9 = dd9[k]; da = dda[k]; db = ddb[k]; dc = ddc[k]; 
	    dd = ddd[k]; de = dde[k];

	    f0 = ff0[k]; f1 = ff1[k]; f2 = ff2[k]; f3 = ff3[k]; f4 = ff4[k]; f5 = ff5[k]; 
	    f6 = ff6[k]; f7 = ff7[k]; f8 = ff8[k]; f9 = ff9[k]; fa = ffa[k]; fb = ffb[k]; 
	    fc = ffc[k]; fd = ffd[k]; fe = ffe[k]; 
	 
	    nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k];

	    dn_dx = (n[d1] - n[d2])/2;
	    dn_dy = (n[d3] - n[d4])/2;
     	    if (mask[k] == 26) {
       		if (subMask[k] == 0) dn_dz= phi11;
       		if (subMask[k] == 1) dn_dz= phi12;
      		if (subMask[k] == 2) dn_dz= phi180; 
		del2n = n[d1] + n[d2] + n[d3] + n[d4] - 4*nn + 2*(n[d5] - dn_dz - nn);     
     	    }
     	    else if (mask[k] == 27) {
       		if (subMask[k] == 0) dn_dz= - phi11; // check the minus sign
       		if (subMask[k] == 1) dn_dz= - phi12; // check the minus sign
      		if (subMask[k] == 2) dn_dz= -phi180; 
      		del2n = n[d1] + n[d2] + n[d3] + n[d4] - 4*nn + 2*(n[d6] + dn_dz - nn);     
     	    }
	    else {
		dn_dz = (n[d5] - n[d6])/2;
	        del2n = (n[d1] + n[d2] + n[d3] + n[d4] + n[d5] + n[d6] - 6*nn);
	    }

	    equilibrium();
	    
	    tau = tauGas + (tauLiquid - tauGas)*(1+tanh(3*(nn-nc)))/2;

	    fn0[k] = f0 + (fe0 - f0)/tau;
	    fn1[k] = f1 + (fe1 - f1)/tau + n[k]*GG[1];
	    fn2[k] = f2 + (fe2 - f2)/tau + n[k]*GG[2];
	    fn3[k] = f3 + (fe3 - f3)/tau + n[k]*GG[3];
	    fn4[k] = f4 + (fe4 - f4)/tau + n[k]*GG[4];
	    fn5[k] = f5 + (fe5 - f5)/tau + n[k]*GG[5];
	    fn6[k] = f6 + (fe6 - f6)/tau + n[k]*GG[6];
	    fn7[k] = f7 + (fe7 - f7)/tau + n[k]*GG[7];
	    fn8[k] = f8 + (fe8 - f8)/tau + n[k]*GG[8];
	    fn9[k] = f9 + (fe9 - f9)/tau + n[k]*GG[9];
	    fna[k] = fa + (fea - fa)/tau + n[k]*GG[10];
	    fnb[k] = fb + (feb - fb)/tau + n[k]*GG[11];
	    fnc[k] = fc + (fec - fc)/tau + n[k]*GG[12];
	    fnd[k] = fd + (fed - fd)/tau + n[k]*GG[13];
	    fne[k] = fe + (fee - fe)/tau + n[k]*GG[14];

	    /*fn0[k] = f0 + (fe0 - f0)/tau;
	    fn1[k] = f1 + (fe1 - f1)/tau + n[k]*GG[1] * (1+tanh(3*(nn-nc)))/2;
	    fn2[k] = f2 + (fe2 - f2)/tau + n[k]*GG[2] * (1+tanh(3*(nn-nc)))/2;
	    fn3[k] = f3 + (fe3 - f3)/tau + n[k]*GG[3] * (1+tanh(3*(nn-nc)))/2;
	    fn4[k] = f4 + (fe4 - f4)/tau + n[k]*GG[4] * (1+tanh(3*(nn-nc)))/2;
	    fn5[k] = f5 + (fe5 - f5)/tau + n[k]*GG[5] * (1+tanh(3*(nn-nc)))/2;
	    fn6[k] = f6 + (fe6 - f6)/tau + n[k]*GG[6] * (1+tanh(3*(nn-nc)))/2;
	    fn7[k] = f7 + (fe7 - f7)/tau + n[k]*GG[7] * (1+tanh(3*(nn-nc)))/2;
	    fn8[k] = f8 + (fe8 - f8)/tau + n[k]*GG[8] * (1+tanh(3*(nn-nc)))/2;
	    fn9[k] = f9 + (fe9 - f9)/tau + n[k]*GG[9] * (1+tanh(3*(nn-nc)))/2;
	    fna[k] = fa + (fea - fa)/tau + n[k]*GG[10] * (1+tanh(3*(nn-nc)))/2;
	    fnb[k] = fb + (feb - fb)/tau + n[k]*GG[11] * (1+tanh(3*(nn-nc)))/2;
	    fnc[k] = fc + (fec - fc)/tau + n[k]*GG[12] * (1+tanh(3*(nn-nc)))/2;
	    fnd[k] = fd + (fed - fd)/tau + n[k]*GG[13] * (1+tanh(3*(nn-nc)))/2;
	    fne[k] = fe + (fee - fe)/tau + n[k]*GG[14] * (1+tanh(3*(nn-nc)))/2;*/
	    
    }
  }
}

// ===============================================================
/* This function performs the streaming step of LB simulation   */
// ===============================================================

void propagation()
{

  for( k = 0 ; k < N ; k++) {
    if (mask[k] != 28) {
      ff0[k] = fn0[k]; 
      ff1[k] = fn1[k]; 
      ff2[k] = fn2[k]; 
      ff3[k] = fn3[k]; 
      ff4[k] = fn4[k]; 
      ff5[k] = fn5[k]; 
      ff6[k] = fn6[k]; 
      ff7[k] = fn7[k]; 
      ff8[k] = fn8[k]; 
      ff9[k] = fn9[k]; 
      ffa[k] = fna[k]; 
      ffb[k] = fnb[k]; 
      ffc[k] = fnc[k]; 
      ffd[k] = fnd[k]; 
      ffe[k] = fne[k]; 
    }
  }
  for( k = 0 ; k < N ; k++) {
    if (mask[k] != 28) {
      if (dd1[k] != -1) ff1[dd1[k]] = fn1[k]; 
      if (dd2[k] != -1) ff2[dd2[k]] = fn2[k]; 
      if (dd3[k] != -1) ff3[dd3[k]] = fn3[k]; 
      if (dd4[k] != -1) ff4[dd4[k]] = fn4[k]; 
      if (dd5[k] != -1) ff5[dd5[k]] = fn5[k]; 
      if (dd6[k] != -1) ff6[dd6[k]] = fn6[k]; 
      if (dd7[k] != -1) ff7[dd7[k]] = fn7[k]; 
      if (dd8[k] != -1) ff8[dd8[k]] = fn8[k]; 
      if (dd9[k] != -1) ff9[dd9[k]] = fn9[k]; 
      if (dda[k] != -1) ffa[dda[k]] = fna[k]; 
      if (ddb[k] != -1) ffb[ddb[k]] = fnb[k]; 
      if (ddc[k] != -1) ffc[ddc[k]] = fnc[k]; 
      if (ddd[k] != -1) ffd[ddd[k]] = fnd[k]; 
      if (dde[k] != -1) ffe[dde[k]] = fne[k]; 
    }
  }

}

// ============================================================
/* Boundary condition used to complete partial distribution  */
/* function. Only applies for 2-D flow. Need to change for   */
/* 3-D flow.                                                 */ 
// ============================================================

void applyBoundaryConditions()
{
  
  double r = 0;

  for(k = k1 ; k < k2 ; k++) {

    switch(mask[k]) {
      
    case 0:
      break;

    case 26:
      r=drand48();
#ifdef NORANDOM
      r=0;
#endif
      if (r < 0.25) {
	ff7[k] = ff8[k];
	ff5[k] = ff6[k];
	ffb[k] = 0.5*(ff1[k]-ff2[k]+ff3[k]-ff4[k])+ffc[k];
	ff9[k] = 0.5*(-ff3[k]+ff4[k])+ffa[k];
	ffd[k] = 0.5*(-ff1[k]+ff2[k])+ffe[k];
      }
      else if ((r >= 0.25) && (r < 0.5)) {
	ffb[k] = 0.5*(ff1[k]-ff2[k])+ffc[k];
	ffd[k] = 0.5*(ff3[k]-ff4[k]-ff1[k]+ff2[k])+ffe[k];
	ff9[k] = ffa[k];
	ff7[k] = 0.5*(-ff3[k]+ff4[k])+ff8[k];
	ff5[k] = ff6[k];
      }
      else if ((r >= 0.5) && (r < 0.75)) {
	ffb[k] = ffc[k];
	ff7[k] = 0.5*(-ff1[k]+ff2[k]-ff3[k]+ff4[k])+ff8[k];
	ff9[k] = 0.5*(ff1[k]-ff2[k])+ffa[k];
	ffd[k] = 0.5*(ff3[k]-ff4[k])+ffe[k];
	ff5[k] = ff6[k];
      }
      else if (r >= 0.75) {
	ffd[k] = ffe[k]; 
	ff7[k] = 0.5*(-ff1[k]+ff2[k])+ff8[k];
	ffb[k] = 0.5*(ff3[k]-ff4[k])+ffc[k]; 
	ff9[k] = 0.5*(ff1[k]-ff2[k]-ff3[k]+ff4[k])+ffa[k];
	ff5[k] = ff6[k];
      } 
      ff0[k]+= -(ff5[k]+ff7[k]+ffb[k]+ff9[k]+ffd[k])+
	        (fn6[k]+fn8[k]+fnc[k]+fna[k]+fne[k]);      
      break;
 
    case 27:
      r=drand48();
#ifdef NORANDOM
      r=0;
#endif
      if (r < 0.25) {
	ffa[k] = 0.5*(ff3[k]-ff4[k])+ff9[k];
	ff6[k] = ff5[k];
	ff8[k] = ff7[k];
	ffe[k] = 0.5*(ff1[k]-ff2[k])+ffd[k];
	//ffc[k] = -0.5*(-ff3[k]+ff4[k]-ff1[k]+ff2[k])+ffb[k]; 
	ffc[k] = 0.5*(-ff3[k]+ff4[k]-ff1[k]+ff2[k])+ffb[k]; 
      }
      else if ((r >= 0.25) && (r < 0.5)) {
	ff6[k] = ff5[k];
	ffc[k] = 0.5*(-ff1[k]+ff2[k])+ffb[k];
	ff8[k] = 0.5*(ff3[k]-ff4[k])+ff7[k];
	ffe[k] = 0.5*(ff1[k]-ff2[k]-ff3[k]+ff4[k])+ffd[k]; 
	ffa[k] = ff9[k];
      }
      else if ((r >= 0.5) && (r < 0.75)) {
	ff6[k] = ff5[k];
	ffc[k] = ffb[k];
	ffe[k] = 0.5*(-ff3[k]+ff4[k])+ffd[k];
	ffa[k] = 0.5*(-ff1[k]+ff2[k])+ff9[k];
	ff8[k] = 0.5*(ff1[k]-ff2[k]+ff3[k]-ff4[k])+ff7[k];
      }
      else if (r >= 0.75) {
	ff6[k] = ff5[k];
	ffe[k] = ffd[k]; 
	ffc[k] = 0.5*(-ff3[k]+ff4[k])+ffb[k];
	ff8[k] = 0.5*(ff1[k]-ff2[k])+ff7[k]; 
	ffa[k] = ff9[k]+0.5*(-ff1[k]+ff2[k]+ff3[k]-ff4[k]);
      }
      ff0[k]+= -(ff6[k]+ff8[k]+ffc[k]+ffe[k]+ffa[k])+
	        (fn5[k]+fn7[k]+fn9[k]+fnb[k]+fnd[k]); 
      break;

    default:
      std::cout << "A case is not defined (" << (int) (mask[k]) << ")." << std::endl;
      break;
    }
  }
}

// =======================================================================
/* This function synchronises the partial density fn* at the y-boundary */
// =======================================================================

void synchronisefn()
{

  int xLattice = 0;
  int tempk = 0;

  for (k = 0; k < N; k++) {
	
    computeCoordinate();
    xLattice = int(k/((float)(LZ*LY)));

    if (yk == 0) {

      tempk = zk + 2*LZ + xLattice*LZ*LY; 	  

      fn0[k] = fn0[tempk];
      fn1[k] = fn1[tempk];
      fn2[k] = fn2[tempk];
      fn3[k] = fn4[tempk];
      fn4[k] = fn3[tempk];
      fn5[k] = fn5[tempk];
      fn6[k] = fn6[tempk];
      fn7[k] = fnd[tempk];
      fn8[k] = fne[tempk];
      fn9[k] = fnb[tempk];
      fna[k] = fnc[tempk];
      fnb[k] = fn9[tempk];
      fnc[k] = fna[tempk];
      fnd[k] = fn7[tempk];
      fne[k] = fn8[tempk];

    }

    if (yk == LY - 1) {

      tempk = zk + (LY-3)*LZ + xLattice*LZ*LY; 	  

      fn0[k] = fn0[tempk];
      fn1[k] = fn1[tempk];
      fn2[k] = fn2[tempk];
      fn3[k] = fn4[tempk];
      fn4[k] = fn3[tempk];
      fn5[k] = fn5[tempk];
      fn6[k] = fn6[tempk];
      fn7[k] = fnd[tempk];
      fn8[k] = fne[tempk];
      fn9[k] = fnb[tempk];
      fna[k] = fnc[tempk];
      fnb[k] = fn9[tempk];
      fnc[k] = fna[tempk];
      fnd[k] = fn7[tempk];
      fne[k] = fn8[tempk];

    }
	
  }

}

// =======================================================================
/* This function synchronises the partial density ff* at the y-boundary */
// =======================================================================

void synchroniseff()
{

  int xLattice = 0;
  int tempk = 0;

  for (k = 0; k < N; k++) {
	
    computeCoordinate();
    xLattice = int(k/((float)(LZ*LY)));

    if (yk == 0) {

      tempk = zk + 2*LZ + xLattice*LZ*LY; 	  

      ff0[k] = ff0[tempk];
      ff1[k] = ff1[tempk];
      ff2[k] = ff2[tempk];
      ff3[k] = ff4[tempk];
      ff4[k] = ff3[tempk];
      ff5[k] = ff5[tempk];
      ff6[k] = ff6[tempk];
      ff7[k] = ffd[tempk];
      ff8[k] = ffe[tempk];
      ff9[k] = ffb[tempk];
      ffa[k] = ffc[tempk];
      ffb[k] = ff9[tempk];
      ffc[k] = ffa[tempk];
      ffd[k] = ff7[tempk];
      ffe[k] = ff8[tempk];

    }

    if (yk == LY - 1) {

      tempk = zk + (LY-3)*LZ + xLattice*LZ*LY; 	  

      ff0[k] = ff0[tempk];
      ff1[k] = ff1[tempk];
      ff2[k] = ff2[tempk];
      ff3[k] = ff4[tempk];
      ff4[k] = ff3[tempk];
      ff5[k] = ff5[tempk];
      ff6[k] = ff6[tempk];
      ff7[k] = ffd[tempk];
      ff8[k] = ffe[tempk];
      ff9[k] = ffb[tempk];
      ffa[k] = ffc[tempk];
      ffb[k] = ff9[tempk];
      ffc[k] = ffa[tempk];
      ffd[k] = ff7[tempk];
      ffe[k] = ff8[tempk];

    }
	
  }

}


// ============================================================
/* This function computes the 3-D coordinate of index k      */
/* Filling order: zk, yk, xk                                 */
// ============================================================

void computeCoordinate()
{  

  xk=int(k/(float) (LZ*LY));
  yk=int((k-xk*LZ*LY)/(float) LZ);
  zk=k-xk*LZ*LY-yk*LZ;

#ifdef PARALLEL
  xk+=myPE*LX/nbPE-1;
#endif

}

// ==============================================================
/* This function computes the velocities of each lattice point */
// ==============================================================

void computeMomenta()
{
  for( k = k1 ; k < k2 ; k++) {
    f0 = ff0[k]; f1 = ff1[k]; f2 = ff2[k]; f3 = ff3[k]; f4 = ff4[k]; 
    f5 = ff5[k]; f6 = ff6[k]; f7 = ff7[k]; f8 = ff8[k]; f9 = ff9[k]; 
    fa = ffa[k]; fb = ffb[k]; fc = ffc[k]; fd = ffd[k]; fe = ffe[k]; 
    
    nn = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + fa + fb + fc + fd + fe;
    n[k] = nn;
    uxs[k] = (f1 + f7 + fa + fc + fd - f2 - f8 - f9 - fb - fe)/nn; 
    uys[k] = (f3 + f7 + f9 + fc + fe - f4 - f8 - fa - fb - fd)/nn;
    uzs[k] = (f5 + f7 + f9 + fb + fd - f6 - f8 - fa - fc - fe)/nn;
  }
}

// =============================================================
/* This function computes the total free energy of the system */
/* There is a mistake in the gradient term - corrected Halim  */
/* Also, density gradient must be calculated for each k       */
// =============================================================

void computeFreeEnergy()
{
   energy = 0.0;
   freeE = 0.0;
   interfaceE = 0.0;
   surfaceE = 0.0;
   
   for( k = k1 ; k < k2 ; k++) {

// ###########################################################
/* Bulk term                                                */
// ###########################################################
     
     nue= n[k]/nc - 1.0;
     freeE += pc*(nue+1)*(nue+1)*(nue*nue-2*nue+3-2*beta*gam);
     
// ###########################################################
/* Interface term                                           */
// ###########################################################
     
     dn_dx = ( n[dd1[k]] - n[dd2[k]] ) / 2;
     dn_dy = ( n[dd3[k]] - n[dd4[k]] ) / 2;
     if (mask[k] == 26) {
       if (subMask[k] == 0) dn_dz= phi11;
       if (subMask[k] == 1) dn_dz= phi12;
       if (subMask[k] == 2) dn_dz= phi180; 
     }
     else if (mask[k] == 27) {
       if (subMask[k] == 0) dn_dz= - phi11; // check the minus sign
       if (subMask[k] == 1) dn_dz= - phi12; // check the minus sign
       if (subMask[k] == 2) dn_dz= - phi180; 
     }
     else {
       dn_dz = ( n[dd5[k]] - n[dd6[k]] ) / 2;
     }
     
     interfaceE += kappa/2.0*(dn_dx*dn_dx + dn_dy*dn_dy + dn_dz*dn_dz);

// ###########################################################
/* Surface term                                             */
// ###########################################################

     if (mask[k] == 26 || mask[k] == 27) {
       if (subMask[k] == 0) surfaceE += phi11*kappa*n[k]; 
       else if (subMask[k] == 1) surfaceE += phi12*kappa*n[k];                 
       else surfaceE += phi180*kappa*n[k];
     }
   
#ifdef SYMMETRY
     
     computeCoordinate();

     if (yk > 2 || yk < LY -3) {

       freeE += pc*(nue+1)*(nue+1)*(nue*nue-2*nue+3-2*beta*gam);
       interfaceE += kappa/2.0*(dn_dx*dn_dx + dn_dy*dn_dy + dn_dz*dn_dz);
       if (mask[k] == 26 || mask[k] == 27) {
	 if (subMask[k] == 0) surfaceE += phi11*kappa*n[k]; 
	 else if (subMask[k] == 1) surfaceE += phi12*kappa*n[k];                 
	 else surfaceE += phi180*kappa*n[k];
       }

     }


#endif
  
   }

   energy = freeE + interfaceE + surfaceE;
   
#ifdef PARALLEL

  double reducedEnergy = 0.0;
  MPI_Reduce(&energy,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  energy=reducedEnergy;

  reducedEnergy = 0.0;
  MPI_Reduce(&freeE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  freeE=reducedEnergy;

  reducedEnergy = 0.0;
  MPI_Reduce(&interfaceE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  interfaceE=reducedEnergy;

  reducedEnergy = 0.0;
  MPI_Reduce(&surfaceE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  surfaceE=reducedEnergy;

#endif

}

// ============================================================
/* This function returns the total density of the system     */
// ============================================================

double getTotalDensity()
{
  double result=0.0;
  for(k = k1 ; k < k2 ; k++) {

    result+=ff0[k]+ff1[k]+ff2[k]+ff3[k]+ff4[k]+ff5[k]+ff6[k]+ff7[k]+ff8[k]+ff9[k]+ffa[k]+ffb[k]+ffc[k]+ffd[k]+ffe[k]; 

#ifdef SYMMETRY
    computeCoordinate();
    if (yk > 2 || yk < LY -3)
      result += ff0[k]+ff1[k]+ff2[k]+ff3[k]+ff4[k]+ff5[k]+ff6[k]+ff7[k]+ff8[k]+ff9[k]+ffa[k]+ffb[k]+ffc[k]+ffd[k]+ffe[k]; 
#endif

  }

#ifdef PARALLEL
  double reducedResult=0.0;
  MPI_Reduce(&result,&reducedResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return reducedResult;
#else
  return result;
#endif

}

// ============================================================
/* This function returns the droplet's total density         */
// ============================================================

double getDropletDensity()
{
  double result = 0.0;
  double tempn = 0.0;
  for(k = k1 ; k < k2 ; k++) {
    tempn = ff0[k]+ff1[k]+ff2[k]+ff3[k]+ff4[k]+ff5[k]+ff6[k]+ff7[k]+ff8[k]+ff9[k]+ffa[k]+ffb[k]+ffc[k]+ffd[k]+ffe[k]; 
    if ( tempn > ((nFluid + nGas)/2) )
      result+= tempn;

#ifdef SYMMETRY
    computeCoordinate();
    if (yk > 2 || yk < LY -3) { 
      if ( tempn > ((nFluid + nGas)/2) )
	result+= tempn;
    }
#endif

  }

#ifdef PARALLEL
  double reducedResult=0.0;
  MPI_Reduce(&result,&reducedResult,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  return reducedResult;
#else
  return result;
#endif

}

// =============================================================
/* This function computes the centre of mass of the system    */
/* It can be used in both 2-D and 3-D                         */
/* It uses comX,comY,comZ to save the centre of mass position */
/* and comVX, comVY, and comVZ for the velocity               */
// =============================================================

void computeCentreOfMass()
{

// ##########################################################
/* Define temporary variables                              */
/* 1) leftboundary = 1 and rightboundary = 1 if droplet is */
/* positioned across the periodic boundary position        */
/* 2) overlap = LX is added to the position of fluid       */
/* lattice points near the left boundary                   */
// ##########################################################

  int leftboundary = 0, rightboundary = 0;
  int overlap = 0;
  double comXtemp = comX  - ( ((int) (int)comX / LX) ) * LX;
  double dropletDensity = 0;
	
// ##########################################################
/* Reinitialise the variables in use                       */
// ##########################################################

  comX = 0; comY = 0; comZ = 0;
  comVX = 0; comVY = 0; comVZ = 0;
  computeMomenta();

// ##########################################################
/* Scan x = 0 and x = LX - 1 to find out if droplet cross  */
/* the periodic boundary condition                         */
// ##########################################################

#ifdef PARALLEL
  if (myPE == 0) {
    i = 1;
    for (j = 0; j < LY; j++) {
      for (h = 0; h < LZ; h++) {
	k = h + j*LZ + i*LY*LZ;	  
	if ( n[k] > ((nFluid + nGas)/2) ) leftboundary = 1;
      }
    }
  }
  
  if (myPE == nbPE - 1) {
    i = i2 - 2;
    for (j = 0; j < LY; j++) {
      for (h = 0; h < LZ; h++) {
	k = h + j*LZ + i*LY*LZ;	  
	if ( n[k] > ((nFluid + nGas)/2) ) rightboundary = 1;
      }
    }
  }
  
  MPI_Bcast(&leftboundary, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rightboundary, 1, MPI_INT, nbPE-1, MPI_COMM_WORLD);
	
#else	
  i = 0;
  for (j = 0; j < LY; j++) {
    for (h = 0; h < LZ; h++) {
      k = h + j*LZ + i*LY*LZ;	  
      if ( n[k] > ((nFluid + nGas)/2) ) leftboundary = 1;
    }	
  }

  i = LX - 1;	
  for (j = 0; j < LY; j++) {
    for (h = 0; h < LZ; h++) {              	
      k = h + j*LZ + i*LY*LZ;	  
      if ( n[k] > ((nFluid + nGas)/2) ) rightboundary = 1;
    }
  }

#endif

// ##########################################################
/* Compute the COM variables                               */
// ##########################################################

  for (k = k1; k < k2; k++) {
    computeCoordinate();
    if ( n[k] > ((nFluid + nGas)/2) ) {
      dropletDensity += n[k];
      if ( leftboundary == 1 && rightboundary == 1 && xk < (int)(LX/2) ) overlap = LX;
      else overlap = 0;
      comX += n[k]*(xk+overlap); comY += n[k]*yk; comZ += n[k]*zk;
      comVX += n[k]*uxs[k]; comVY += n[k]*uys[k]; comVZ += n[k]*uzs[k];

#ifdef SYMMETRY
      if (yk > 2 || yk < LY -3) { 
	dropletDensity += n[k];
	comX += n[k]*(xk+overlap); comY += n[k]*(2-yk); comZ += n[k]*zk;
	comVX += n[k]*uxs[k]; comVY -= n[k]*uys[k]; comVZ += n[k]*uzs[k];
      }
#endif

    }
  }
	
#ifdef PARALLEL

  double temporary = 0.0;
  MPI_Reduce(&dropletDensity,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  dropletDensity=temporary;
  
  temporary = 0.0;
  MPI_Reduce(&comX,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comX=temporary;
  
  temporary = 0.0;
  MPI_Reduce(&comY,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comY=temporary;

  temporary = 0.0;
  MPI_Reduce(&comZ,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comZ=temporary;

  temporary = 0.0;
  MPI_Reduce(&comVX,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comVX=temporary;

  temporary = 0.0;
  MPI_Reduce(&comVY,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comVY=temporary;

  temporary = 0.0;
  MPI_Reduce(&comVZ,&temporary,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  comVZ=temporary;
  
#endif

  if (myPE == 0) {
  
    comX = comX/dropletDensity; comY = comY/dropletDensity; comZ = comZ/dropletDensity;
    comVX = comVX/dropletDensity; comVY = comVY/dropletDensity; comVZ = comVZ/dropletDensity;
    if (comX > (double) LX) comX -= (double) LX;
	
// ##########################################################
/* Housekeeping: How many times have the droplet cross     */
/* the periodic boundary condition? Add crossing*LX to     */
/* comX to ensure comX is monotonically increasing.        */                                     
// ##########################################################

    if (comXtemp - comX > dropletR) crossing++;
    comX += crossing * LX;

// ##########################################################
/* Save the COM variables                                  */
/* check std::ios::app                                          */
/* save in comX.m comY.m comZ.m comVX.m comVY.m comVZ.m    */
// ##########################################################

    std::ofstream fout1("comX.m", std::ios::app);
    fout1 << comX << std::endl;
    fout1.close();
	
    std::ofstream fout2("comY.m", std::ios::app);
    fout2 << comY << std::endl;
    fout2.close();
	
    std::ofstream fout3("comZ.m", std::ios::app);
    fout3 << comZ << std::endl;
    fout3.close();
	
    std::ofstream fout4("comVX.m", std::ios::app);
    fout4 << comVX << std::endl;
    fout4.close();
	
    std::ofstream fout5("comVY.m", std::ios::app);
    fout5 << comVY << std::endl;
    fout5.close();

    std::ofstream fout6("comVZ.m", std::ios::app);
    fout6 << comVZ << std::endl;
    fout6.close();
	
  }
  
}

// ============================================================
/* This function writes the density of the lattice points    */ 
/* into matlab m-file: dt(:,:,h+1)= [...]                    */ 
/* filename: ddt.m where t = number of iterations            */
/* and h is zk (i.e. density is written in x-y plane slices) */
// ============================================================

void writeDensityFile()
{

#ifdef PARALLEL
  MPI_Gather(&n[k1],k2-k1,MPI_DOUBLE,nGlobal,k2-k1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#else
  nGlobal=n;
#endif

  if (myPE == 0) {
    String fileName("dd");
    fileName.concat((int) t);
    fileName.concat(".m");
    std::ofstream file(fileName.get());
    file.precision(4);
    for( h = 0 ; h < LZ ; h++) {   
      file << "d" << t << "(:,:," << h+1 << ")=[" << std::endl;
      for( i = 0 ; i < LX ; i++) {
	for( j = 0 ; j < LY ; j++) {
	  k = h + j*LZ + i*LY*LZ;	  
	  file << nGlobal[k] << " ";
	}
	file << std::endl;
      }
      file << "];" << std::endl;
    }
    file.close();
  }

}

// ============================================================
/* This function writes the density of the lattice points at */
/* certain xz plane (y = iy = j) into dyt=[...]              */
/* filename: dyjtt.m. Second t is the number of iteration    */
// ============================================================

void writeYPlanDensityFile(const int iy)
{
  j=iy;

#ifdef PARALLEL
  MPI_Gather(&n[k1],k2-k1,MPI_DOUBLE,nGlobal,k2-k1,MPI_DOUBLE,
	     0,MPI_COMM_WORLD);
#else
  nGlobal=n;
#endif

  if (myPE == 0) {
    String fileName("dy");
    fileName.concat((int) j);
    fileName.concat("t");
    fileName.concat((int) t);
    fileName.concat(".m");
    std::ofstream file(fileName.get());
    file.precision(4);
    file << "dy" << t << "=[" << std::endl;
    for( i = 0 ; i < LX ; i++) {
      for( h = 0 ; h < LZ ; h++) {   
	k = h + j*LZ + i*LY*LZ;	
	file << nGlobal[k] << " ";
      }
      file << std::endl;
    }
    file << "];" << std::endl;
    file.close();
  }
}

// ============================================================
/* This function writes the density of the lattice points at */
/* the substrate boundary (z = 0) into st.m : st=[...]       */
// ============================================================

void writeDensityFileOnSubstrate()
{

#ifdef PARALLEL
  MPI_Gather(&n[k1],k2-k1,MPI_DOUBLE,nGlobal,k2-k1,MPI_DOUBLE,
	     0,MPI_COMM_WORLD);
#else
  nGlobal=n;
#endif

  if (myPE == 0) {
    String fileName("s");
    fileName.concat((int) t);
    fileName.concat(".m");
    std::ofstream file(fileName.get());
    file.precision(4);
    h=0;
    file << "s" << t << "=[" << std::endl;
    for( i = 0 ; i < LX ; i++) {
      for( j = 0 ; j < LY ; j++) {
	k = h + j*LZ + i*LY*LZ;	  
	file << nGlobal[k] << " ";
      }
      file << std::endl;
    }
    file << "];" << std::endl;
    file.close();
  }
}

// ============================================================
/* This function writes the velocities of the lattice points */
/* at certain xz plane (y = j = iy) into uxYjtt.m: uxt=[...] */
/* for x-direction and uzYjtt.m: uzt=[...] for z-direction   */
// ============================================================

void writeYPlanXVelocityFile(const int iy)
{
  j=iy;

#ifdef PARALLEL
  MPI_Gather(&uxs[k1],k2-k1,MPI_DOUBLE,nGlobal,k2-k1,MPI_DOUBLE,
	     0,MPI_COMM_WORLD);
#else
  nGlobal=uxs;
#endif

  if (myPE == 0) {
    String fileName("uxY");
    fileName.concat((int) j);
    fileName.concat("t");
    fileName.concat((int) t);
    fileName.concat(".m");
    std::ofstream file(fileName.get());
    file.precision(4);
    file << "ux" << t << "=[" << std::endl;
    for( i = 0 ; i < LX ; i++) {
      for( h = 0 ; h < LZ ; h++) {   
	k = h + j*LZ + i*LY*LZ;	
	file << nGlobal[k] << " ";
      }
      file << std::endl;
    }
    file << "];" << std::endl;
    file.close();
  }


#ifdef PARALLEL
  MPI_Gather(&uzs[k1],k2-k1,MPI_DOUBLE,nGlobal,k2-k1,MPI_DOUBLE,
	     0,MPI_COMM_WORLD);
#else
  nGlobal=uzs;
#endif

  if (myPE == 0) {
    String fileName("uzY");
    fileName.concat((int) j);
    fileName.concat("t");
    fileName.concat((int) t);
    fileName.concat(".m");
    std::ofstream file(fileName.get());
    file.precision(4);
    file << "uz" << t << "=[" << std::endl;
    for( i = 0 ; i < LX ; i++) {
      for( h = 0 ; h < LZ ; h++) {   
	k = h + j*LZ + i*LY*LZ;	
	file << nGlobal[k] << " ";
      }
      file << std::endl;
    }
    file << "];" << std::endl;
    file.close();
  }

}


#ifdef PARALLEL
// ============================================================
void exchangeBoundaries()
{
  // Sends to neighbors
  //
  // Left
  //
  // The required fields are: 2,8,9,11(b),14(e).
  MPI_Bsend(&(fn2[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,0,MPI_COMM_WORLD);
  MPI_Bsend(&(fn8[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,1,MPI_COMM_WORLD);
  MPI_Bsend(&(fn9[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,2,MPI_COMM_WORLD);
  MPI_Bsend(&(fnb[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,3,MPI_COMM_WORLD);
  MPI_Bsend(&(fne[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,4,MPI_COMM_WORLD);
  //MPI_Bsend(&(fn2[k1]),1,leftFieldsType,leftNeighbor,0,MPI_COMM_WORLD);

  //
  // Right
  //
  // The required fields are: 1,7,10(a),12(c),13(d).
  MPI_Bsend(&(fn1[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,0,MPI_COMM_WORLD);
  MPI_Bsend(&(fn7[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,1,MPI_COMM_WORLD);
  MPI_Bsend(&(fna[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,2,MPI_COMM_WORLD);
  MPI_Bsend(&(fnc[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,3,MPI_COMM_WORLD);
  MPI_Bsend(&(fnd[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,4,MPI_COMM_WORLD);
  //MPI_Bsend(&(fn1[k2-LY*LZ]),1,rightFieldsType,rightNeighbor,0,MPI_COMM_WORLD);
 
  // Receives from neighbors
  //
  // Right
  //
  MPI_Recv(&(fn2[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,0,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fn8[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,1,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fn9[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,2,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fnb[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,3,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fne[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,4,MPI_COMM_WORLD,&status);
  //MPI_Recv(&(fn2[k2]),1,leftFieldsType,rightNeighbor,0,MPI_COMM_WORLD,&status);

  //
  // Left
  //
  MPI_Recv(&(fn1[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,0,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fn7[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,1,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fna[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,2,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fnc[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,3,MPI_COMM_WORLD,&status);
  MPI_Recv(&(fnd[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,4,MPI_COMM_WORLD,&status);
  //MPI_Recv(&(fn1[0]),1,rightFieldsType,leftNeighbor,0,MPI_COMM_WORLD,&status);

}

// ============================================================
void exchangeDensitiesAtBoundaries()
{

  // Sends
  //
  MPI_Bsend(&(n[k1]),LZ*LY,MPI_DOUBLE,leftNeighbor,0,MPI_COMM_WORLD);
  //
  MPI_Bsend(&(n[k2-LY*LZ]),LZ*LY,MPI_DOUBLE,rightNeighbor,0,MPI_COMM_WORLD);

  // Receives
  // 
  MPI_Recv(&(n[k2]),LZ*LY,MPI_DOUBLE,rightNeighbor,0,MPI_COMM_WORLD,&status);
  //
  MPI_Recv(&(n[0]),LZ*LY,MPI_DOUBLE,leftNeighbor,0,MPI_COMM_WORLD,&status);
  //

}
#endif

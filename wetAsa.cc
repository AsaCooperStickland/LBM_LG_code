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
  
  // ASA: This needs chamging 
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

    /*double a = (xk-dropletCenterX)*(xk-dropletCenterX)  + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-dropletCenterZ)*(zk-dropletCenterZ);
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
/* Allocate memory for variables (ASA: there are now 27, previously 15)                          */
// ##########################################################

  ff0 = new double[N]; 
  
  ffa1 = new double[N]; ffa2= new double[N]; ffa3 = new double[N]; 
  ffa4 = new double[N]; ffa5= new double[N]; ffa6 = new double[N]; 
  
  ffb1 = new double[N]; ffb2 = new double[N]; ffb3 = new double[N]; 
  ffb4 = new double[N]; ffb5 = new double[N]; ffb6 = new double[N]; 
  ffb7 = new double[N]; ffb8 = new double[N]; ffb9 = new double[N];
  ffba = new double[N]; ffbb = new double[N]; ffbc = new double[N];

  ffc1 = new double[N]; ffc2= new double[N]; ffc3 = new double[N]; 
  ffc4 = new double[N]; ffc5= new double[N]; ffc6 = new double[N];
  ffc7= new double[N]; ffc8 = new double[N];


  fn0 = new double[N]; 
  
  fna1 = new double[N]; fna2= new double[N]; fna3 = new double[N]; 
  fna4 = new double[N]; fna5= new double[N]; fna6 = new double[N]; 
  
  fnb1 = new double[N]; fnb2 = new double[N]; fnb3 = new double[N]; 
  fnb4 = new double[N]; fnb5 = new double[N]; fnb6 = new double[N]; 
  fnb7 = new double[N]; fnb8 = new double[N]; fnb9 = new double[N];
  fnba = new double[N]; fnbb = new double[N]; fnbc = new double[N];

  fnc1 = new double[N]; fnc2= new double[N]; fnc3 = new double[N]; 
  fnc4 = new double[N]; fnc5= new double[N]; fnc6 = new double[N];
  fnc7= new double[N]; fnc8 = new double[N];

  n = new double[N];

  uxs = new double[N];
  uys = new double[N];
  uzs = new double[N];

  dd0 = new long[N]; 
  
  dda1 = new long[N]; dda2= new long[N]; dda3 = new long[N]; 
  dda4 = new long[N]; dda5= new long[N]; dda6 = new long[N]; 
  
  ddb1 = new long[N]; ddb2 = new long[N]; ddb3 = new long[N]; 
  ddb4 = new long[N]; ddb5 = new long[N]; ddb6 = new long[N]; 
  ddb7 = new long[N]; ddb8 = new long[N]; ddb9 = new long[N];
  ddba = new long[N]; ddbb = new long[N]; ddbc = new long[N];

  ddc1 = new long[N]; ddc2= new long[N]; ddc3 = new long[N]; 
  ddc4 = new long[N]; ddc5= new long[N]; ddc6 = new long[N];
  ddc7= new long[N]; ddc8 = new long[N]; 

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

    // ##########################################################
    /* ASA: The bit below is still for 15 directions, needs changing!                       */
    // ##########################################################
    
    dda1[k] = zk + yk*LZ + r*LZ*LY; 	  
    dda2[k] = zk + yk*LZ + l*LZ*LY; 
    dda3[k] = zk + u*LZ + xkLattice*LZ*LY; 
    dda4[k] = zk + d*LZ + xkLattice*LZ*LY; 
    dda5[k] = w + yk*LZ + xkLattice*LZ*LY; 
    dda6[k] = q + yk*LZ + xkLattice*LZ*LY; 
    ddb1[k] = zk + u*LZ + r*LZ*LY; 
    ddb2[k] = zk + d*LZ + r*LZ*LY; 
    ddb3[k] = zk + u*LZ + l*LZ*LY; 
    ddb4[k] = zk + d*LZ + l*LZ*LY; 
    ddb5[k] = w + yk*LZ + r*LZ*LY; 
    ddb6[k] = q + yk*LZ + r*LZ*LY; 
    ddb7[k] = w + yk*LZ + l*LZ*LY; 
    ddb8[k] = q + yk*LZ + l*LZ*LY; 
    ddb9[k] = w + u*LZ + xkLattice*LZ*LY;
    ddba[k] = q + u*LZ + xkLattice*LZ*LY;
    ddbb[k] = w + d*LZ + xkLattice*LZ*LY;
    ddbc[k] = q + d*LZ + xkLattice*LZ*LY;
    ddc1[k] = w + u*LZ + r*LZ*LY;
    ddc2[k] = q + u*LZ + r*LZ*LY;
    ddc3[k] = w + d*LZ + r*LZ*LY;
    ddc4[k] = q + d*LZ + r*LZ*LY;
    ddc5[k] = w + u*LZ + l*LZ*LY;
    ddc6[k] = q + u*LZ + l*LZ*LY;
    ddc7[k] = w + d*LZ + l*LZ*LY;
    ddc8[k] = q + d*LZ + l*LZ*LY;

//ASA: Need to figure out this bit
    if (zk == 0)
      {
	dda6[k] = -1;
	ddb6[k] = -1;
	ddb8[k] = -1;
  ddba[k] = -1;
	ddbc[k] = -1;
	ddc2[k] = -1;
  ddc4[k] = -1;
  ddc6[k] = -1;
  ddc8[k] = -1;
      }
    
    if (zk == LZ - 1)
      {
	dda5[k] = -1;
	ddb5[k] = -1;
	ddb7[k] = -1;
  ddb9[k] = -1;
	ddbb[k] = -1;
	ddc1[k] = -1;
  ddc3[k] = -1;
  ddc5[k] = -1;
  ddc7[k] = -1;
      } 

#ifdef PARALLEL
    if (xkLattice == 0)
    {
	dda2[k] = -1;
	ddb1[k] = -1;
	ddb2[k] = -1;
	ddb7[k] = -1;
	ddb8[k] = -1;
  ddc5[k] = -1;
  ddc6[k] = -1;
  ddc7[k] = -1;
  ddc8[k] = -1;	
    }    
    if (xkLattice == i2 - 1)
    {
	dda1[k] = -1;
	ddb1[k] = -1;
	ddb2[k] = -1;
	ddb5[k] = -1;
	ddb6[k] = -1;
  ddc1[k] = -1;
  ddc2[k] = -1;
  ddc3[k] = -1;
  ddc4[k] = -1;    
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
    ffa1[k] = fea1; ffa2[k] = fea2; ffa3[k] = fea3; ffa4[k] = fea4; ffa5[k] = fea5; ffa6[k] = fea6; 
    ffb1[k] = feb1; ffb2[k] = feb2; ffb3[k] = feb3; ffb4[k] = feb4; ffb5[k] = feb5; ffb6[k] = feb6; 
    ffb7[k] = feb7; ffb8[k] = feb8; ffb9[k] = feb9; ffba[k] = feba; ffbb[k] = febb; ffbc[k] = febc;
    ffc1[k] = fec1; ffc2[k] = fec2; ffc3[k] = fec3; ffc4[k] = fec4; ffc5[k] = fec5; ffc6[k] = fec6; 
    ffc7[k] = fec7; ffc8[k] = fec8;
  
  }

// ###########################################################
/* Velocity vector definition                               */
// ###########################################################

// This MAYBE needs changing
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
  uxy = ux*uy; uxz = ux*uz; uyz = uy*uz; 
  u2=ux2+uy2+uz2;
  double nue= nn/nc - 1.0;

  // Should this be a cont in the hh file? 
  cs2 = 0.333333333;
  cs2_div = 1.0/(cs2);
  cs4_div = 1.0/(cs2*cs2);
  cs6_div = 1.0/(cs2*cs2*cs2);
  
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
  // ============================================================
  /* ASA: have completely changed what follows, see original code for what it was. Corresponds to Eq. 6 in the paper 
  Need to allocate memory to these!!       */
  // ============================================================
  fea1 = w1*(nn*(1.0 + cs2_div*ux + ux2*(0.5*cs2_div + 0.16667*cs6_div*ux) - u2*0.5*(cs2_div + cs4_div*ux)));
  fea2 = w1*(nn*(1.0 - cs2_div*ux + ux2*(0.5*cs2_div - 0.16667*cs6_div*ux) - u2*0.5*(cs2_div - cs4_div*ux)));
  fea3 = w1*(nn*(1.0 + cs2_div*uy + uy2*(0.5*cs2_div + 0.16667*cs6_div*uy) - u2*0.5*(cs2_div + cs4_div*uy)));
  fea4 = w1*(nn*(1.0 - cs2_div*uy + uy2*(0.5*cs2_div - 0.16667*cs6_div*uy) - u2*0.5*(cs2_div - cs4_div*uy)));
  fea5 = w1*(nn*(1.0 + cs2_div*uz + uz2*(0.5*cs2_div + 0.16667*cs6_div*uz) - u2*0.5*(cs2_div + cs4_div*uz)));
  fea6 = w1*(nn*(1.0 - cs2_div*uz + uz2*(0.5*cs2_div - 0.16667*cs6_div*uz) - u2*0.5*(cs2_div - cs4_div*uz)));

  feb1 = w2*(nn*(1 + cs2_div*(ux + uy) + (ux2 + 2*uxy + uy2 )*(0.5*cs2_div + 0.16667*cs6_div*(ux + uy)) 
  - u2*0.5*(cs2_div + cs4_div*(ux + uy)))); 
  feb2 = w2*(nn*(1 + cs2_div*(ux - uy) + (ux2 - 2*uxy + uy2 )*(0.5*cs2_div + 0.16667*cs6_div*(ux - uy)) 
  - u2*0.5*(cs2_div + cs4_div*(ux - uy)))); 
  feb3 = w2*(nn*(1 + cs2_div*(-ux + uy) + (ux2 - 2*uxy + uy2 )*(0.5*cs2_div + 0.16667*cs6_div*(-ux + uy)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux + uy)))); 
  feb4 = w2*(nn*(1 + cs2_div*(-ux - uy) + (ux2 + 2*uxy + uy2 )*(0.5*cs2_div + 0.16667*cs6_div*(-ux - uy)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux - uy)))); 
  feb5 = w2*(nn*(1 + cs2_div*(ux + uz) + (ux2 + 2*uxz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(ux + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux + uz)))); 
  feb6 = w2*(nn*(1 + cs2_div*(ux - uz) + (ux2 - 2*uxz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(ux - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux - uz)))); 
  feb7 = w2*(nn*(1 + cs2_div*(-ux + uz) + (ux2 - 2*uxz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(-ux + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux + uz)))); 
  feb8 = w2*(nn*(1 + cs2_div*(-ux - uz) + (ux2 + 2*uxz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(-ux - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux - uz)))); 
  feb9 = w2*(nn*(1 + cs2_div*(uy + uz) + (uy2 + 2*uyz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(uy + uz)))); 
  feba = w2*(nn*(1 + cs2_div*(uy - uz) + (uy2 - 2*uyz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(uy - uz)))); 
  febb = w2*(nn*(1 + cs2_div*(-uy + uz) + (uy2 - 2*uyz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(-uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-uy + uz)))); 
  febc = w2*(nn*(1 + cs2_div*(-uy - uz) + (uy2 + 2*uyz + uz2 )*(0.5*cs2_div + 0.16667*cs6_div*(-uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-uy - uz))));
 
  fec1 = w3*(nn*(1 + cs2_div*(ux + uy + uz) + (u2 + 2*(uxy + uxz + uyz))*(0.5*cs2_div + 0.16667*cs6_div*(ux + uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux + uy + uz)))); 
  fec2 = w3*(nn*(1 + cs2_div*(ux + uy - uz) + (u2 + 2*(uxy - uxz - uyz))*(0.5*cs2_div + 0.16667*cs6_div*(ux + uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux + uy - uz)))); 
  fec3 = w3*(nn*(1 + cs2_div*(ux - uy + uz) + (u2 + 2*(-uxy + uxz - uyz))*(0.5*cs2_div + 0.16667*cs6_div*(ux - uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux - uy + uz)))); 
  fec4 = w3*(nn*(1 + cs2_div*(ux - uy - uz) + (u2 + 2*(-uxy - uxz + uyz))*(0.5*cs2_div + 0.16667*cs6_div*(ux - uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(ux - uy - uz)))); 
  fec5 = w3*(nn*(1 + cs2_div*(-ux + uy + uz) + (u2 + 2*(-uxy - uxz + uyz))*(0.5*cs2_div + 0.16667*cs6_div*(-ux + uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux + uy + uz)))); 
  fec6 = w3*(nn*(1 + cs2_div*(-ux + uy - uz) + (u2 + 2*(-uxy + uxz - uyz))*(0.5*cs2_div + 0.16667*cs6_div*(-ux + uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux + uy - uz)))); 
  fec7 = w3*(nn*(1 + cs2_div*(-ux - uy + uz) + (u2 + 2*(uxy - uxz - uyz))*(0.5*cs2_div + 0.16667*cs6_div*(-ux - uy + uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux - uy + uz)))); 
  fec8 = w3*(nn*(1 + cs2_div*(-ux - uy - uz) + (u2 + 2*(uxy + uxz + uyz))*(0.5*cs2_div + 0.16667*cs6_div*(-ux - uy - uz)) 
  - u2*0.5*(cs2_div + cs4_div*(-ux - uy - uz))));

  // ASA: Need to change this as well
  fe0 = nn - fea1 - fea2 - fea3 - fea4 - 
        fea5 - fea6 - feb1 - feb2 - feb3 - feb4 - feb5 - feb6 - feb7 - feb8 - 
        feb9 - feba - febb - febc - fec1 - fec2 - fec3 - fec4 - fec5 - fec6 - fec7 - fec8;
      
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


	    // ASA: is this needed other than for x,y,z i.e. a1,a2...,a6? 
	    d1 = dda1[k]; d2 = dda2[k]; d3 = dda3[k]; d4 = dda4[k]; d5 = dda5[k]; d6 = dda6[k];
	   
        // ASA: I've added these, NEED TO CREATE THE VARIABLES IN HH FILE still. 
        // Some I don't htink are needed (not used in dxx, dxy etc., will delete at some point'
		
		d11 = dda1[dda1[k]]; d13 = dda1[dda3[k]]; d14 = dda1[dda4[k]]; d15 = dda1[dda5[k]]; d16 = dda1[dda6[k]];   
		d21 = dda2[dda1[k]]; d22 = dda2[dda2[k]]; d23 = dda2[dda3[k]]; d24 = dda2[dda4[k]]; d25 = dda2[dda5[k]]; d26 = dda2[dda6[k]]; 
		d33 = dda3[dda3[k]]; d35 = dda3[dda5[k]]; d36 = dda3[dda6[k]];  
		d43 = dda4[dda3[k]]; d44 = dda4[dda4[k]]; d45 = dda4[dda5[k]]; d46 = dda4[dda6[k]]; 
		d55 = dda5[dda5[k]]; 
		d66 = dda6[dda6[k]];   
		
		f0 = ff0[k]; fa1 = ffa1[k]; fa2 = ffa2[k]; fa3 = ffa3[k]; fa4 = ffa4[k]; fa5 = ffa5[k]; fa6 = ffa6[k]; 
		fb1 = ffb1[k]; fb2 = ffb2[k]; fb3 = ffb3[k]; fb4 = ffb4[k]; fb5 = ffb5[k]; fb6 = ffb6[k]; fb7 = ffb7[k]; 
		fb8 = ffb8[k]; fb9 = ffb9[k]; fba = ffba[k]; fbb = ffbb[k]; fbc = ffbc[k]; 
		fc1 = ffc1[k]; fc2 = ffc2[k]; fc3 = ffc3[k]; fc4 = ffc4[k]; fc5 = ffc5[k]; fc6 = ffc6[k]; fc7 = ffc7[k]; fc8 = ffc8[k]; 
	 
	    nn = n[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k];

        // ASA: I've added the 2nd order partials here. WHY is dn_dz not always defined??

	    dn_dx = (n[d1] - n[d2])/2;
	    dn_dy = (n[d3] - n[d4])/2;
		
		dxx = (n[d1] + n[d2] - 2 * nn);
		dxy = (n[d13] - n[d14] - n[d23] + n[d24]) / 4;
		dxz = (n[d15] - n[d16] - n[d25] + n[d26]) / 4;
		dyy = (n[d3] + n[d4] - 2 * nn);
		dyz = (n[d35] - n[d36] - n[d45] + n[d46]) / 4;
		dzz = (n[d5] + n[d6] - 2 * nn);

		// ASA: This is dxxx + dyyyy + dzzz

		dthird = (-n[d11] + 2 * n[d1] - 2 * n[d2] + n[d22] - n[d33] + 2 * n[d3] - 2 * n[d4] + n[d44] - n[d55] + 2 * n[d5] - 2 * n[d6] + n[d66]) / 2;
        // Should dn_dz go here??
     	    if (mask[k] == 26) {
       		if (subMask[k] == 0) dn_dz= phi11;
       		if (subMask[k] == 1) dn_dz= phi12;
      		if (subMask[k] == 2) dn_dz= phi180; 

 

      // Need to create these variabls in the hh file!!! apr, bpr, rg are from the EoS in the paper. T is tempearture
		apr = 9/49;
		bpr = 2/21;
		rg = 1;
		T2 = T*T;
		T3 = T2*T;
		T4 = T3*T;
		T5 = T4*T;
		T6 = T5*T;
		// is it faster to define nn2 = nn*nn, nn3 = nn*nn2 etc???
		// Might be an issue with placement of this stuff
		dpr_dn = rg*T + (rg*bpr*(2.416231612436E+05*T - 2.57967008222E+07*T2 + 1.14238248744E+09*T3 - 2.68526135606E+10*T4 + 3.533712827532E+11*T5
		- 2.468605110492E+12*T6) - 2 * apr)*nn + rg*bpr*(-1.2404676319134E+05*T + 1.3260848313711E+07*T2 - 5.878424929332E+08*T3 + 1.3833808099563E+10*T4
		- 1.8230179071144E+11*T5 + 1.2756315594687E+12*T6)*nn*nn + rg*bpr*(4.1456868307976E+02*T - 5.9411760240720E+04*T2 + 3.4719996817435E+06*T3
		- 1.0533330468360E+08*T4 + 1.7490689935796E+09*T5 - 1.5090810131432E+10*T6)*nn*nn*nn + rg*bpr*(9.1264500505120E+03*T - 9.6910045234700E+05*T2
		+ 4.2609615223305E+07*T3 - 9.9308940538450E+08*T4 + 1.2941593445460E+10*T5 - 8.9421093774750E+10*T6)*nn*nn*nn*nn + rg*bpr*(-1.0858736022620E+03*T
		+ 1.1535804904722E+05*T2 - 5.0760202853244E+06*T3 + 1.1843360670144E+08*T4 - 1.5455391750738E+09*T5 + 1.0697254852476E+10*T6)*nn*nn*nn*nn*nn;
		ka = 0.0; // ka is the k in the paper
      
     
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

		dfe0 = fe0 - f0; dfea1 = fea1 - fa1; dfea2 = fea2 - fa2; dfea3 = fea3 - fa3; dfea4 = fea4 - fa4; dfea5 = fea5 - fa5; dfea6 = fea6 - fa6;
		dfeb1 = feb1 - fb1; dfeb2 = feb2 - fb2; dfeb3 = feb3 - fb3; dfeb4 = feb4 - fb4; dfeb5 = feb5 - fb5; dfeb6 = feb6 - fb6;
		dfeb7 = feb7 - fb7; dfeb8 = feb8 - fb8; dfeb9 = feb9 - fb9; dfeba = feba - fba; dfebb = febb - fbb; dfebc = febc - fbc;
		dfec1 = fec1 - fc1; dfec2 = fec2 - fc2; dfec3 = fec3 - fc3; dfec4 = fec4 - fc4; dfec5 = fec5 - fc5; dfec6 = fec6 - fc6; dfec7 = fec7 - fc7; dfec8 = fec8 - fc8;

		// ASA: I added these temporary fe variables
		tfe0 = fe0; tfea1 = fea1; tfea2 = fea2; tfea3 = fea3; tfea4 = fea4; tfea5 = fea5; tfea6 = fea6;
		tfeb1 = feb1; tfeb2 = feb2; tfeb3 = feb3; tfeb4 = feb4; tfeb5 = feb5; tfeb6 = feb6;
		tfeb7 = feb7; tfeb8 = feb8; tfeb9 = feb9; tfeba = feba; tfebb = febb; tfebc = febc;
		tfec1 = fec1; tfec2 = fec2; tfec3 = fec3; tfec4 = fec4; tfec5 = fec5; tfec6 = fec6; tfec7 = fec7; tfec8 = fec8;

		lw = log(w);
		lw1 = log(w1);
		lw2 = log(w2);
		lw3 = log(w3);

		al0 = 2.0;

		for (h = 0; h < 10; h++) {
			entropy = (f0 + al0*dfe0)*log(f0 + al0*dfe0) - f0*log(f0) - al0*dfe0*lw + (fa1 + al0*dfea1)*log(fa1 + al0*dfea1) - fa1*log(fa1) - al0*dfea1*lw1 +
			(fa2 + al0*dfea2)*log(fa2 + al0*dfea2) - fa2*log(fa2) - al0*dfea2*lw1 + (fa3 + al0*dfea3)*log(fa3 + al0*dfea3) - fa3*log(fa3) - al0*dfea3*lw1 +
			(fa4 + al0*dfea4)*log(fa4 + al0*dfea4) - fa4*log(fa4) - al0*dfea4*lw1 + (fa5 + al0*dfea5)*log(fa5 + al0*dfea5) - fa5*log(fa5) - al0*dfea5*lw1 +
			(fa6 + al0*dfea6)*log(fa6 + al0*dfea6) - fa6*log(fa6) - al0*dfea6*lw1 + (fb1 + al0*dfeb1)*log(fb1 + al0*dfeb1) - fb1*log(fb1) - al0*dfeb1*lw2 +
			(fb2 + al0*dfeb2)*log(fb2 + al0*dfeb2) - fb2*log(fb2) - al0*dfeb2*lw2 + (fb3 + al0*dfeb3)*log(fb3 + al0*dfeb3) - fb3*log(fb3) - al0*dfeb3*lw2 +
			(fb4 + al0*dfeb4)*log(fb4 + al0*dfeb4) - fb4*log(fb4) - al0*dfeb4*lw2 + (fb5 + al0*dfeb5)*log(fb5 + al0*dfeb5) - fb5*log(fb5) - al0*dfeb5*lw2 +
			(fb6 + al0*dfeb6)*log(fb6 + al0*dfeb6) - fb6*log(fb6) - al0*dfeb6*lw2 + (fb7 + al0*dfeb7)*log(fb7 + al0*dfeb7) - fb7*log(fb7) - al0*dfeb7*lw2 +
			(fb8 + al0*dfeb8)*log(fb8 + al0*dfeb8) - fb8*log(fb8) - al0*dfeb8*lw2 + (fb9 + al0*dfeb9)*log(fb9 + al0*dfeb9) - fb9*log(fb9) - al0*dfeb9*lw2 +
			(fba + al0*dfeba)*log(fba + al0*dfeba) - fba*log(fba) - al0*dfeba*lw2 + (fbb + al0*dfebb)*log(fbb + al0*dfebb) - fbb*log(fbb) - al0*dfebb*lw2 +
			(fbc + al0*dfebc)*log(fbc + al0*dfebc) - fbc*log(fbc) - al0*dfebc*lw2 + (fc1 + al0*dfec1)*log(fc1 + al0*dfec1) - fc1*log(fc1) - al0*dfec1*lw3 +
			(fc2 + al0*dfec2)*log(fc2 + al0*dfec2) - fc2*log(fc2) - al0*dfec2*lw3 + (fc3 + al0*dfec3)*log(fc3 + al0*dfec3) - fc3*log(fc3) - al0*dfec3*lw3 +
			(fc4 + al0*dfec4)*log(fc4 + al0*dfec4) - fc4*log(fc4) - al0*dfec5*lw3 + (fc5 + al0*dfec5)*log(fc5 + al0*dfec5) - fc5*log(fc5) - al0*dfec5*lw3 +
			(fc6 + al0*dfec6)*log(fc6 + al0*dfec6) - fc6*log(fc6) - al0*dfec6*lw3 + (fc7 + al0*dfec7)*log(fc7 + al0*dfec7) - fc7*log(fc7) - al0*dfec7*lw3 +
			(fc8 + al0*dfec8)*log(fc8 + al0*dfec8) - fc8*log(fc8) - al0*dfec8*lw3;

			dentropy = (al0*dfe0)*log(f0 + al0*dfe0) + (fa1 + al0*dfea1)*log(fa1 + al0*dfea1) +
			(al0*dfea2)*log(fa2 + al0*dfea2) + (fa3 + al0*dfea3)*log(fa3 + al0*dfea3) +
			(al0*dfea4)*log(fa4 + al0*dfea4) + (fa5 + al0*dfea5)*log(fa5 + al0*dfea5) +
			(al0*dfea6)*log(fa6 + al0*dfea6) + (fb1 + al0*dfeb1)*log(fb1 + al0*dfeb1) +
			(al0*dfeb2)*log(fb2 + al0*dfeb2) + (fb3 + al0*dfeb3)*log(fb3 + al0*dfeb3) +
			(al0*dfeb4)*log(fb4 + al0*dfeb4) + (fb5 + al0*dfeb5)*log(fb5 + al0*dfeb5) +
			(al0*dfeb6)*log(fb6 + al0*dfeb6) + (fb7 + al0*dfeb7)*log(fb7 + al0*dfeb7) +
			(al0*dfeb8)*log(fb8 + al0*dfeb8) + (fb9 + al0*dfeb9)*log(fb9 + al0*dfeb9) +
			(al0*dfeba)*log(fba + al0*dfeba) + (fbb + al0*dfebb)*log(fbb + al0*dfebb) +
			(al0*dfebc)*log(fbc + al0*dfebc) + (fc1 + al0*dfec1)*log(fc1 + al0*dfec1) +
			(al0*dfec2)*log(fc2 + al0*dfec2) + (fc3 + al0*dfec3)*log(fc3 + al0*dfec3) +
			(al0*dfec4)*log(fc4 + al0*dfec4) + (fc5 + al0*dfec5)*log(fc5 + al0*dfec5) +
			(al0*dfec6)*log(fc6 + al0*dfec6) + (fc7 + al0*dfec7)*log(fc7 + al0*dfec7) +
			(al0*dfec8)*log(fc8 + al0*dfec8) + dfe0 + dfea1 + dfea2 + dfea3 + dfea4 +
			dfea5 + dfea6 + dfeb1 + dfeb2 + dfeb3 + dfeb4 + dfeb5 + dfeb6 + dfeb7 + dfeb8 +
			dfeb9 + dfeba + dfebb + dfebc + dfec1 + dfec2 + dfec3 + dfec4 + dfec5 + dfec6 + dfec7 + dfec8;

			al1 = al0 - entropy / dentropy;
			al0 = al1;

		}
		// Here I add force term, derivations are in labbook. have assumed dt = 1.0
		ux = ux + (cs2*dn_dx - dn_dx*dpr_dn + ka*(dn_dx*(-dxx + del2n + dxy + dxz) + dn_dy*(dxx + dxy) + dn_dz*(dxx + dxz) + nn*dthird)) / nn;
		uy = uy + (cs2*dn_dy - dn_dy*dpr_dn + ka*(dn_dy*(-dyy + del2n + dxy + dyz) + dn_dx*(dyy + dxy) + dn_dz*(dyy + dyz) + nn*dthird)) / nn;
		uz = uz + (cs2*dn_dz - dn_dz*dpr_dn + ka*(dn_dz*(-dzz + del2n + dxz + dyz) + dn_dx*(dzz + dxz) + dn_dy*(dzz + dyz) + nn*dthird)) / nn;
		// Set new equibm after changing the speed 
		equilibrium();

		tau = tauGas + (tauLiquid - tauGas)*(1 + tanh(3 * (nn - nc))) / 2;

		alm = 1 - al1 / tau;

		fn0[k] = alm*(f0 - tfe0) + fe0;

	    fna1[k] = alm*(fa1 - tfea1) + fea1;
	    fna2[k] = alm*(fa2 - tfea2) + fea2;
	    fna3[k] = alm*(fa3 - tfea3) + fea3;
		fna4[k] = alm*(fa4 - tfea4) + fea4;
	    fna5[k] = alm*(fa5 - tfea5) + fea5;
	    fna6[k] = alm*(fa6 - tfea6) + fea6;

	    fnb1[k] = alm*(fb1 - tfeb1) + feb1;
	    fnb2[k] = alm*(fb2 - tfeb2) + feb2;
	    fnb3[k] = alm*(fb3 - tfeb3) + feb3;
		fnb4[k] = alm*(fb4 - tfeb4) + feb4;
	    fnb5[k] = alm*(fb5 - tfeb5) + feb5;
	    fnb6[k] = alm*(fb6 - tfeb6) + feb6;
		fnb7[k] = alm*(fb7 - tfeb7) + feb7;
	    fnb8[k] = alm*(fb8 - tfeb8) + feb8;
	    fnb9[k] = alm*(fb9 - tfeb9) + feb9;
		fnba[k] = alm*(fba - tfeba) + feba;
	    fnbb[k] = alm*(fbb - tfebb) + febb;
	    fnbc[k] = alm*(fbc - tfebc) + febc;
		
		fnc1[k] = alm*(fc1 - tfec1) + fec1;
	    fnc2[k] = alm*(fc2 - tfec2) + fec2;
	    fnc3[k] = alm*(fc3 - tfec3) + fec3;
		fnc4[k] = alm*(fc4 - tfec4) + fec4;
	    fnc5[k] = alm*(fc5 - tfec5) + fec5;
	    fnc6[k] = alm*(fc6 - tfec6) + fec6;
		fnc7[k] = alm*(fc7 - tfec7) + fec7;
	    fnc8[k] = alm*(fc8 - tfec8) + fec8;

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
      ffa1[k] = fna1[k]; 
      ffa2[k] = fna2[k]; 
      ffa3[k] = fna3[k]; 
      ffa4[k] = fna4[k]; 
      ffa5[k] = fna5[k]; 
      ffa6[k] = fna6[k]; 
  
      ffb1[k] = fnb1[k]; 
      ffb2[k] = fnb2[k]; 
      ffb3[k] = fnb3[k]; 
      ffb4[k] = fnb4[k]; 
      ffb5[k] = fnb5[k]; 
      ffb6[k] = fnb6[k]; 
      ffb7[k] = fnb7[k]; 
      ffb8[k] = fnb8[k]; 
      ffb9[k] = fnb9[k]; 
      ffba[k] = fnba[k]; 
      ffbb[k] = fnbb[k]; 
      ffbc[k] = fnbc[k]; 

      ffc1[k] = fnc1[k]; 
      ffc2[k] = fnc2[k]; 
      ffc3[k] = fnc3[k]; 
      ffc4[k] = fnc4[k]; 
      ffc5[k] = fnc5[k]; 
      ffc6[k] = fnc6[k]; 
      ffc7[k] = fnc7[k]; 
      ffc8[k] = fnc8[k];

    }
  }
  for( k = 0 ; k < N ; k++) {
    if (mask[k] != 28) {
      if (dda1[k] != -1) ffa1[dda1[k]] = fna1[k]; 
      if (dda2[k] != -1) ffa2[dda2[k]] = fna2[k]; 
      if (dda3[k] != -1) ffa3[dda3[k]] = fna3[k]; 
      if (dda4[k] != -1) ffa4[dda4[k]] = fna4[k]; 
      if (dda5[k] != -1) ffa5[dda5[k]] = fna5[k]; 
      if (dda6[k] != -1) ffa6[dda6[k]] = fna6[k];

      if (ddb1[k] != -1) ffb1[ddb1[k]] = fnb1[k]; 
      if (ddb2[k] != -1) ffb2[ddb2[k]] = fnb2[k]; 
      if (ddb3[k] != -1) ffb3[ddb3[k]] = fnb3[k]; 
      if (ddb4[k] != -1) ffb4[ddb4[k]] = fnb4[k]; 
      if (ddb5[k] != -1) ffb5[ddb5[k]] = fnb5[k]; 
      if (ddb6[k] != -1) ffb6[ddb6[k]] = fnb6[k]; 
      if (ddb7[k] != -1) ffb7[ddb7[k]] = fnb7[k]; 
      if (ddb8[k] != -1) ffb8[ddb8[k]] = fnb8[k]; 
      if (ddb9[k] != -1) ffb9[ddb9[k]] = fnb9[k]; 
      if (ddba[k] != -1) ffba[ddba[k]] = fnba[k]; 
      if (ddbb[k] != -1) ffbb[ddbb[k]] = fnbb[k]; 
      if (ddbc[k] != -1) ffbc[ddbc[k]] = fnbc[k]; 

      if (ddc1[k] != -1) ffc1[ddc1[k]] = fnc1[k]; 
      if (ddc2[k] != -1) ffc2[ddc2[k]] = fnc2[k]; 
      if (ddc3[k] != -1) ffc3[ddc3[k]] = fnc3[k]; 
      if (ddc4[k] != -1) ffc4[ddc4[k]] = fnc4[k]; 
      if (ddc5[k] != -1) ffc5[ddc5[k]] = fnc5[k]; 
      if (ddc6[k] != -1) ffc6[ddc6[k]] = fnc6[k]; 
      if (ddc7[k] != -1) ffc7[ddc7[k]] = fnc7[k]; 
      if (ddc8[k] != -1) ffc8[ddc8[k]] = fnc8[k];
      
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

    // This should work in 27 directions now
    f0 = ff0[k]; fa1 = ffa1[k]; fa2 = ffa2[k]; fa3 = ffa3[k]; fa4 = ffa4[k]; fa5 = ffa5[k]; fa6 = ffa6[k]; 
    fb1 = ffb1[k]; fb2 = ffb2[k]; fb3 = ffb3[k]; fb4 = ffb4[k]; fb5 = ffb5[k]; fb6 = ffb6[k]; fb7 = ffb7[k]; 
    fb8 = ffb8[k]; fb9 = ffb9[k]; fba = ffba[k]; fbb = ffbb[k]; fbc = ffbc[k]; 
    fc1 = ffc1[k]; fc2 = ffc2[k]; fc3 = ffc3[k]; fc4 = ffc4[k]; fc5 = ffc5[k]; fc6 = ffc6[k]; fc7 = ffc7[k]; fc8 = ffc8[k];  
    
    nn = f0 + fa1 + fa2 + fa3 + fa4 + fa5 + fa6 + fb1 + fb2 + fb3 + fb4 + fb5 + fb6 + fb7 + fb8 + fb9 + fba + fbb + fbc +
    fc1 + fc2 + fc3 + fc4 + fc5 + fc6 + fc7 + fc8;
    n[k] = nn;
    uxs[k] = (fa1 + fb1 + fb2 + fb5 + fb6 + fc1 + fc2 + fc3 + fc4  - fa2 - fb3 - fb4 - fb7 - fb8 - fc5 - fc6 - fc7 - fc8)/nn; 
    uys[k] = (fa3 + fb1 + fb3 + fb9 + fba + fc1 + fc2 + fc5 + fc6 - fa4 - fb2 - fb4 - fba - fbc - fc3 - fc4 - fc7 - fc8)/nn;
    uzs[k] = (fa5 + fb5 + fb7 + fb9 + fbb + fc1 + fc3 + fc5 + fc7 - fa6 - fb6 - fb8 - fba - fbc - fc2 - fc4 - fc6 - fc8)/nn;
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

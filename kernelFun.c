#include <sys/types.h>
#include <sys/stat.h> // create directory
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernelFun.h"
#include "AnisoFunctions.h"


/* ---- Polynomail kernel  ---- */

void PolyFun (vec3* fieldpos, vec3* sourcepos, double *K) {
  
  vec3 diff;
  double r;
	
  // Compute 1/r
  diff.x = sourcepos->x - fieldpos->x;
  diff.y = sourcepos->y - fieldpos->y;
  diff.z = sourcepos->z - fieldpos->z;
  r      = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
	
  int idof, dof2 = 54;
  double axis;
  for (idof=0; idof < dof2; idof++) {
    if (idof%3 == 0)
      axis = diff.x;
    else if (idof%3 == 1)
      axis = diff.y;
    else
      axis = diff.z;
    
    K[idof] = axis * r;
  }
}


/*********************************
 *    LAPLACIAN KERNEL (1/r)     *
 *********************************/

void LaplacianFun (vec3* fieldpos, vec3* sourcepos, double *K) {
  
  vec3 diff;

  diff.x = sourcepos->x - fieldpos->x;
  diff.y = sourcepos->y - fieldpos->y;
  diff.z = sourcepos->z - fieldpos->z;
  
  *K     = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
}


/*********************************
 *  GAUSSIAN KERNEL (exp(-r^2))  *
 *********************************/

void GaussianFun (vec3* fieldpos, vec3* sourcepos, double *K) {

  vec3 diff;
  double r;
	
  diff.x = sourcepos->x - fieldpos->x;
  diff.y = sourcepos->y - fieldpos->y;
  diff.z = sourcepos->z - fieldpos->z;
  r      = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

  *K = exp(-r*r);
}   


/* ---- Laplacianforce kernel  ---- */

void LapforceFun (vec3* fieldpos, vec3* sourcepos, double *K) {
  
  vec3 diff;
  double rinv;
	
  // Compute 1/r
  diff.x = sourcepos->x - fieldpos->x;
  diff.y = sourcepos->y - fieldpos->y;
  diff.z = sourcepos->z - fieldpos->z;
  rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);

  *K =  diff.x*diff.x*diff.x  * rinv*rinv*rinv*rinv*rinv;
  /*
  int idof, dof2 = 54;
  double axis;
  for (idof=0; idof < dof2; idof++) {
    if (idof%3 == 0)
      axis = diff.x;
    else if (idof%3 == 1)
      axis = diff.y;
    else
      axis = diff.z;
    
    K[idof] = axis*axis*axis * rinv*rinv*rinv*rinv*rinv;
  }
  */
}


/* ---- Anisotropic kernel in dislocation dynamics  ---- */

// global anisotropic parameter
static paraAniso* AnisoParameters = 0;

void AnisoFun (vec3* fieldpos, vec3* sourcepos, double *K) {

  assert(AnisoParameters != 0);
  ChaoStress(fieldpos->x, fieldpos->y, fieldpos->z,
	     sourcepos->x, sourcepos->y, sourcepos->z,
	     AnisoParameters->qMax,
	     AnisoParameters->ind1Re, AnisoParameters->ind2Re,
	     AnisoParameters->ind1Im, AnisoParameters->ind2Im,
	     AnisoParameters->ReFqout, AnisoParameters->ImFqout,
	     AnisoParameters->sizeReFq, AnisoParameters->sizeImFq,
	     AnisoParameters->e12, AnisoParameters->e3,
	     AnisoParameters->C633, AnisoParameters->CdGdx,
	     K);
  
}


static void ReadAnisoParameters();
void AnisoInit() {

  assert(AnisoParameters == 0);
  printf("Initializing anisotropic parameters ...\n");
  AnisoParameters = malloc( sizeof(paraAniso) );
  assert(AnisoParameters != NULL);


  // Highly aniso example: need qMax = 20 for this example
  double c11, c12, c44;
  c11=14.88e2;
  c12=12.22e2;
  c44=9.9e2;
  // c11 *= 1e9, c12 *= 1e9, c44 *= 1e9;

  // Isotropic example : qMax = 1 is OK for this example
  /* c11 = 52.1e2;
     c12 = 20.1e2;
     c44 = 16.0e2;
     qMax = 1;
  */

  
  double *c3 = malloc(sizeof(double)*3);
  c3[0] = c11, c3[1] = c12, c3[2] = c44;

  int qMax = 20;
  double elasA = 2*c44/(c11-c12);
  printf("A=%.3f qMax=%d\n", elasA, qMax);

  
  AnisoParameters->qMax = qMax;
  AnisoParameters->elasA = elasA;
  AnisoParameters->c3 = c3;

  
  double *e3 = malloc(sizeof(double)*3);
  e3[0] = 0.0;
  e3[1] = 0.0;
  e3[2] = 1.0;

  real8_complex *e12 = malloc(sizeof(real8_complex)*3);
  e12[0] = 1.0;
  e12[1] = 0 + 1.0*I;
  e12[2] = 0.0;
    
  AnisoParameters->e12 = &e12[0];
  AnisoParameters->e3  = &e3[0];

  
  double C66[6][6], C[3][3][3][3];
  double (*C633)[3][3] = (double (*)[3][3]) malloc( sizeof(double)*54 );
  double (*CdGdx)[6] = (double (*)[6]) malloc( sizeof(double)*18 );

  GetCbcc(AnisoParameters->c3[0], AnisoParameters->c3[1],
	  AnisoParameters->c3[2], C66, C, C633);
   
  AnisoParameters->C633  = C633;
  AnisoParameters->CdGdx = CdGdx;

  ReadAnisoParameters(AnisoParameters);    
     
}


static void PreComputeAniso ();
void ReadAnisoParameters() {

  // create the directory if not exist
  char *dir = "AnisoFiles/";
  struct stat st = {0};
  if (stat(dir, &st) == -1) {
    printf("Create '%s' dirctory.\n", dir);
    mkdir(dir, 0775);
  }
    
  double elasA = AnisoParameters->elasA;  
  char Ind1Re[50], Ind2Re[50], Ind1Im[50], Ind2Im[50];
  char ReFqout[50], ImFqout[50];
	
  sprintf(Ind1Re,"%sInd1ReA%.2f.bin", dir, elasA);
  sprintf(Ind2Re,"%sInd2ReA%.2f.bin", dir, elasA);
  sprintf(Ind1Im,"%sInd1ImA%.2f.bin", dir, elasA);
  sprintf(Ind2Im,"%sInd2ImA%.2f.bin", dir, elasA);
  sprintf(ReFqout,"%sReFqoutA%.2f.bin", dir, elasA);
  sprintf(ImFqout,"%sImFqoutA%.2f.bin", dir, elasA);

  AnisoParameters->f_Ind1Re = Ind1Re;
  AnisoParameters->f_Ind2Re = Ind2Re;
  AnisoParameters->f_Ind1Im = Ind1Im;
  AnisoParameters->f_Ind2Im = Ind2Im;
  AnisoParameters->f_ReFqout = ReFqout;
  AnisoParameters->f_ImFqout = ImFqout;
    
  FILE *Find1Re, *Find2Re, *Find1Im, *Find2Im, *FReFqout, *FImFqout;
  Find1Re = fopen(Ind1Re, "rb");
  Find2Re = fopen(Ind2Re, "rb");
  Find1Im = fopen(Ind1Im, "rb");
  Find2Im = fopen(Ind2Im, "rb");
  FReFqout = fopen(ReFqout, "rb");
  FImFqout = fopen(ImFqout, "rb");

  if (Find1Re == NULL || Find2Re == NULL || Find1Im == NULL || Find2Im
      == NULL || FReFqout == NULL || FImFqout == NULL) {

    printf("Compute anisotropic kernel parameters ...\n");
    PreComputeAniso(AnisoParameters);

    if (FReFqout != NULL)
      fclose(FReFqout);
    if (FImFqout != NULL)
      fclose(FImFqout);
    if (Find1Re != NULL)
      fclose(Find1Re);
    if (Find2Re != NULL)
      fclose(Find2Re);
    if (Find1Im != NULL)
      fclose(Find1Im);
    if (Find2Im != NULL)
      fclose(Find2Im);
      
  } else {

    /*
      fclose(FReFqout);
      fclose(FImFqout);
      fclose(Find1Re);
      fclose(Find2Re);
      fclose(Find1Im);
      fclose(Find2Im);
    */

    //printf("Reading anisotropic kernel parameters...\n");
	 
    // Read 'sizeReFq' and 'sizeImFq' 
    int sizeReFq, sizeImFq, i=0;
    
    //FReFqout = fopen("PreReFqout.bin", "rb");
    //FImFqout = fopen("PreImFqout.bin", "rb");
    i += fread(&sizeReFq, sizeof(int), 1, FReFqout);
    i += fread(&sizeImFq, sizeof(int), 1, FImFqout);

    //printf("sizeReFq:%d, sizeImFq:%d\n", sizeReFq, sizeImFq);
      
    // Read 'ReFqout' and 'ImFqout'
    double *ReFqout = (double *) malloc(sizeReFq*sizeof(double));
    double *ImFqout = (double *) malloc(sizeImFq*sizeof(double));
    i += fread(ReFqout, sizeof(double), sizeReFq, FReFqout);
    i += fread(ImFqout, sizeof(double), sizeImFq, FImFqout);
    fclose(FReFqout);
    fclose(FImFqout);

    // Read 'ind1Re', 'ind2Re', 'ind1Im', and 'ind2Im'
    int *ind1Re = (int *) malloc(sizeReFq*sizeof(int));
    int *ind2Re = (int *) malloc(sizeReFq*sizeof(int));
    int *ind1Im = (int *) malloc(sizeImFq*sizeof(int));
    int *ind2Im = (int *) malloc(sizeImFq*sizeof(int));
      
    //Find1Re = fopen("PreInd1Re.bin", "rb");
    //Find2Re = fopen("PreInd2Re.bin", "rb");
    //Find1Im = fopen("PreInd1Im.bin", "rb");
    //Find2Im = fopen("PreInd2Im.bin", "rb");
    i += fread(ind1Re, sizeof(int), sizeReFq, Find1Re);
    i += fread(ind2Re, sizeof(int), sizeReFq, Find2Re);
    i += fread(ind1Im, sizeof(int), sizeImFq, Find1Im);
    i += fread(ind2Im, sizeof(int), sizeImFq, Find2Im);
      
    fclose(Find1Re);
    fclose(Find2Re);
    fclose(Find1Im);
    fclose(Find2Im);

    if (i != 2*(1+sizeReFq+sizeImFq) + (sizeReFq+sizeImFq))
      printf("fread() error in ReadAnisoParameters().\n");
      
      
    /*
      AnisoParameters = {.sizeReFq = sizeReFq, .sizeImFq = sizeImFq,
      .ReFqout = ReFqout, .ImFqout = ImFqout, .ind1Re =
      ind1Re, .ind2Re = ind2Re, .ind1Im = ind1Im, .ind2Im
      = ind2Im};
    */
      
    AnisoParameters->sizeReFq = sizeReFq;
    AnisoParameters->sizeImFq = sizeImFq;
    AnisoParameters->ind1Re = ind1Re;
    AnisoParameters->ind2Re = ind2Re;
    AnisoParameters->ind1Im = ind1Im;
    AnisoParameters->ind2Im = ind2Im;
    AnisoParameters->ReFqout = ReFqout;
    AnisoParameters->ImFqout = ImFqout;
  }

}




/* 
 * Begin pre-calculations:
 * The following should be located at the beginning of the code and
 * calculated once.
 * The arrays 
 *   ind1Re,ind2Re,
 *   ind1Im,ind2Im,
 *   ReFqout,ImFqout,
 *   sizeReFq, sizeImFq,
 * are filled in the then used to compute the derivative of the Green's function
 */
void PreComputeAniso () {

  /* Pre-calculations need qMax and elastic constants. */

  int numTheta = 150;
  int numPhi = numTheta;
  int numTrapezoidal = 100;
	
  double *phi = (double *)calloc(1, numPhi * sizeof(double));
  double *theta = (double *)calloc(1, numTheta * sizeof(double));
  double *xx = (double *)malloc(numTheta * sizeof(double));
  double *ww = (double *)malloc(numTheta * sizeof(double));

  int i, j;
  for (i = 0; i < numPhi; i++) {
    phi[i] = (double)i / numPhi * 2 * M_PI;
  }

  gqwp(numTheta, xx, ww);
    
  for (i = 0; i < numTheta; i++) {
    theta[i] = acos(xx[i]);
  }
    
  /* full dimensions are dGdxgrid[numTheta][numPhi][3][3][3] */
  double (*dGdxgrid)[numPhi][3][3][3];
  double dGdx[3][3][3];

  dGdxgrid = (double (*)[numPhi][3][3][3]) 
    calloc(1, numTheta * numPhi * 27 * sizeof(double));
    
  double C66[6][6], C[3][3][3][3], C633[6][3][3];
  GetCbcc(AnisoParameters->c3[0], AnisoParameters->c3[1],
	  AnisoParameters->c3[2], C66, C, C633);
	
  /*
   *      Define derivative of the Green's function
   *      Decompose derivative of the Green's function into
   *      a grid of theta and phi angles.
   */
  for (i = 0; i < numTheta; i++) {
    for (j = 0; j < numPhi; j++) {
      int a, b, c;
      double sphericalVecNorm[3];
        
      sphericalVecNorm[0] = sin(theta[i]) * cos(phi[j]);
      sphericalVecNorm[1] = sin(theta[i]) * sin(phi[j]);
      sphericalVecNorm[2] = cos(theta[i]);
        
      GetGreensFuncDerivatives(numTrapezoidal, sphericalVecNorm,
			       C, dGdx);
        
      for (a = 0; a < 3; a++) {
	for (b = 0; b < 3; b++) {
	  for (c = 0; c < 3; c++) {
	    dGdxgrid[i][j][a][b][c] = dGdx[a][b][c];
	  }
	}
      }
    }
  }
    

  /*
   *          Compute expansion coefficients (g^{lm}_{vpg}) of the
   *          derivative of the Green's function.
   *          Compute table of binomial factorial coefficients Q^{lm}(k)
   */
  int qMax = AnisoParameters->qMax;
  real8_complex glm[2*qMax+2][2*qMax+2][3][6];
  double binomialFactCoeffs[2*qMax+2][2*qMax+2][2*qMax+2];

  GetGLMandbinom(qMax, numTheta, numPhi, phi,
		 (double*)dGdxgrid, xx, ww, glm, binomialFactCoeffs);


  /*
   *          Compute function F = Q^{lm}(k) * g^{lm}_{vpg}.
   *          This function is precalculated once and for all at the beginning 
   *          of the simulation. It does not depend on the geometry of the
   *          segments.
   */

  double maxReF, maxImF;
  double *ReFqPtr, *ImFqPtr;
  ReFqPtr = (double *)calloc(1, (qMax+2)*(qMax+1)*18*sizeof(double));
  ImFqPtr = (double *)calloc(1, (qMax+2)*(qMax+1)*18*sizeof(double));
  PreFac(qMax, (double *)binomialFactCoeffs, (real8_complex *)glm,
	 (double *)ReFqPtr, (double *)ImFqPtr, &maxReF, &maxImF);
    
    
  /* 
   *          CREATE LIST 
   *          Create lists for real and imaginary parts of F to flatten F
   *          F[i][j][l]  => ind1list for l
   *                         ind2list for 6*i + j.
   */
  int qMaxProd=(qMax+1)*(qMax+2);
  int maxA=qMaxProd*3*6;
  int sizeReFq, sizeImFq;
	    
  int *ind1Re, *ind2Re;
  double *ReFqout;
	    
  ind1Re = (int *)calloc(maxA,sizeof(int));
  ind2Re = (int *)calloc(maxA,sizeof(int));
  ReFqout = (double *)calloc(maxA,sizeof(double));
	    
  CreateLists(qMax, maxReF, (double *)ReFqPtr,
	      ind1Re, ind2Re, ReFqout, &sizeReFq);
	    
  ind1Re = (int *)realloc(ind1Re,sizeReFq*sizeof(int));
  ind2Re = (int *)realloc(ind2Re,sizeReFq*sizeof(int));
  ReFqout = (double *)realloc(ReFqout,sizeReFq*sizeof(double));
	    
  int *ind1Im, *ind2Im;
  double *ImFqout;
	    
  ind1Im = (int *)calloc(maxA,sizeof(int));
  ind2Im = (int *)calloc(maxA,sizeof(int));
  ImFqout = (double *)calloc(maxA,sizeof(double));
	    
  CreateLists(qMax, maxImF, (double *)ImFqPtr,
	      ind1Im, ind2Im, ImFqout, &sizeImFq);
	    
  ind1Im = (int *)realloc(ind1Im,sizeImFq*sizeof(int));
  ind2Im = (int *)realloc(ind2Im,sizeImFq*sizeof(int));
  ImFqout = (double *)realloc(ImFqout,sizeImFq*sizeof(double));

  AnisoParameters->sizeReFq = sizeReFq;
  AnisoParameters->sizeImFq = sizeImFq;
  AnisoParameters->ind1Re = ind1Re;
  AnisoParameters->ind2Re = ind2Re;
  AnisoParameters->ind1Im = ind1Im;
  AnisoParameters->ind2Im = ind2Im;
  AnisoParameters->ReFqout = ReFqout;
  AnisoParameters->ImFqout = ImFqout;

  FILE *ptr_file;
  ptr_file = fopen(AnisoParameters->f_ReFqout,"wb");
  fwrite(&sizeReFq, sizeof(int), 1, ptr_file);
  fwrite(ReFqout, sizeof(double), sizeReFq, ptr_file);
  fclose(ptr_file);

  ptr_file = fopen(AnisoParameters->f_ImFqout,"wb");
  assert(ptr_file != NULL);
  fwrite(&sizeImFq, sizeof(int), 1, ptr_file);
  fwrite(ImFqout, sizeof(double), sizeImFq, ptr_file);
  fclose(ptr_file);

  ptr_file = fopen(AnisoParameters->f_Ind1Re,"wb");
  fwrite(ind1Re, sizeof(int), sizeReFq, ptr_file);
  fclose(ptr_file);

  ptr_file = fopen(AnisoParameters->f_Ind2Re,"wb");
  fwrite(ind2Re, sizeof(int), sizeReFq, ptr_file);
  fclose(ptr_file);

  ptr_file = fopen(AnisoParameters->f_Ind1Im,"wb");
  fwrite(ind1Im, sizeof(int), sizeImFq, ptr_file);
  fclose(ptr_file);

  ptr_file = fopen(AnisoParameters->f_Ind2Im,"wb");
  fwrite(ind2Im, sizeof(int), sizeImFq, ptr_file);
  fclose(ptr_file);

  free(dGdxgrid);
  free(theta);
  free(phi);
  free(xx);
  free(ww);
  free(ReFqPtr);
  free(ImFqPtr);;
}


// Free anisotropic parameters
void AnisoClean() {

  if (AnisoParameters->c3 != NULL)
    free(AnisoParameters->c3);
  if (AnisoParameters->e3 != NULL)
    free(AnisoParameters->e3);
  if (AnisoParameters->e12 != NULL)
    free(AnisoParameters->e12);
  if (AnisoParameters->C633 != NULL)
    free(AnisoParameters->C633);
  if (AnisoParameters->CdGdx != NULL)
    free(AnisoParameters->CdGdx);
    
  if (AnisoParameters->ind1Re != NULL)
    free(AnisoParameters->ind1Re);
  if (AnisoParameters->ind2Re != NULL)
    free(AnisoParameters->ind2Re);
  if (AnisoParameters->ind1Im != NULL)
    free(AnisoParameters->ind1Im);
  if (AnisoParameters->ind2Im != NULL)
    free(AnisoParameters->ind2Im);
  if (AnisoParameters->ReFqout != NULL)
    free(AnisoParameters->ReFqout);
  if (AnisoParameters->ImFqout != NULL)
    free(AnisoParameters->ImFqout);

  if (AnisoParameters != NULL)
    free(AnisoParameters);
}
 

void Segment2Point(segT *segment, int segmentSize, double *burg, int
		   numGauss, vec3 *source, double *sourceCharge) {

  int i, j;
  double gpoint[numGauss], gweight[numGauss], xi[3*segmentSize];
  gqwp(numGauss, gpoint, gweight); // Gauss points go from 1 to -1
        
  for (i=0; i<segmentSize; i++) {
    xi[3*i+0] = (segment[i].p2.x - segment[i].p1.x)/2;
    xi[3*i+1] = (segment[i].p2.y - segment[i].p1.y)/2;
    xi[3*i+2] = (segment[i].p2.z - segment[i].p1.z)/2;
  }
        
  for (i=0; i<segmentSize; i++) 
    for (j=0;j<numGauss;j++) {
      source[i*numGauss+j].x = (segment[i].p1.x+segment[i].p2.x)/2 + xi[3*i+0] * gpoint[j];
      source[i*numGauss+j].y = (segment[i].p1.y+segment[i].p2.y)/2 + xi[3*i+1] * gpoint[j];
      source[i*numGauss+j].z = (segment[i].p1.z+segment[i].p2.z)/2 + xi[3*i+2] * gpoint[j];
    }
        
  int k, l, count=0;
  for (l=0; l<segmentSize; l++)
    for (k=0; k<numGauss; k++)
      for (i=0; i<3; i++)
	for (j=0; j<3; j++, count++) {
	  //printf("count = %d\n", count);
	  sourceCharge[count] = gweight[k] * burg[3*l+j] * xi[3*l+i];
	}

}


 /* ---- PBC functions ---- */
void DoCorrectionTable(double *ChebyshevWeightSource, double *ChebyshevWeightField, int n, int2 dof, double Len, double alpha, int lpbc, kernel_t kernel, double *Tkz) {

  kfun_t kfun = kernel.kfun;
  double homogen = kernel.homogen;
  
  char kpbcFilename[50];
  sprintf(kpbcFilename, "Kn%dpbc%da%.1f.out", n, lpbc, alpha);
  FILE *kfile = fopen(kpbcFilename, "r");

  // Create correction table
  if (kfile == NULL) {
    printf("kpbc file: %s does not exist. Creating now ...\n", kpbcFilename);
    CreateTableCorrection(kfun, n, dof, alpha, lpbc, Tkz, kpbcFilename);
    kfile = fopen(kpbcFilename, "r");
    assert(kfile != NULL);
  }
  else
    printf("kpbc file exits. Reading now ...\n");
  
  // Read and scale for elastic constants
  int i, j=0;
  double c3Read[3];
  for (i=0; i<3; i++)
    j += fscanf(kfile, "%lf", c3Read+i);
  assert (j == 3);

  int n3 = n*n*n;
  int n3f = n3 * dof.f, n3s = n3 * dof.s;
  int dof2n6 = n3f * n3s;
  double *KPBC = (double *) malloc( dof2n6 * sizeof(double) );
  for (i=0; i<dof2n6; i++) 
    j += fscanf(kfile, "%lf", KPBC+i);
  fclose(kfile);
  assert( j-3 == dof2n6 );
  
  // Compute stress from PBC
  int incr = 1;
  double beta = 0;
  char trans = 'n';
  double scale = pow(1/Len, homogen); // Scale for the box size

  // check if parameters are exactly scaled
  if ( fabs(AnisoParameters->c3[0]/c3Read[0] -
	    AnisoParameters->c3[1]/c3Read[1]) <1e-6 &&   
       fabs(AnisoParameters->c3[2]/c3Read[2] -
	    AnisoParameters->c3[1]/c3Read[1]) <1e-6 )
      
    scale *= AnisoParameters->c3[0]/c3Read[0];
  else
    printf("Error: elastic constants do not match (scale).\n");

  dgemv_(&trans, &n3f, &n3s, &scale, KPBC,
	 &n3f, ChebyshevWeightSource, &incr, &beta,
	 ChebyshevWeightField, &incr);
    
  free(KPBC), KPBC=NULL;

}


static double AdjustBoxSize(double L, double alpha){
  return L * (1+alpha);
}

/*
 * Function: ComputePeriodicKernel
 * ---------------------------------------------------------------------
 * Forms the matrix that describes the interactions of the computational
 * cell with its periodic images up to lpbc shells.
 */
void CreateTableCorrection(kfun_t kfun, int n, int2 dof, double alpha, int lpbc, double *Tkz, char* PBCfilename) {	    

   
  int n3 = n*n*n, i, j, k, l, m, l1, l2, l3, count;              
	
  int dofn3_source = dof.s*n3, dofn3_field = dof.f*n3;
  int dof2n6 = dofn3_source * dofn3_field;
	
  double nodes[n], pi=M_PI;
  vec3 nodepos[n3], cshift, vtmp;

  double L = AdjustBoxSize(1.0, alpha);
  
  // Compute the Chebyshev nodes of T_n(x)
  for (m=0;m<n;m++)
    nodes[m] = cos(pi*((double)m+0.5)/(double)n) * L;

  // Compute the locations of the field points
  count = 0;
  for (l1=0;l1<n;l1++) {
    vtmp.x = 0.5*nodes[l1];
    for (l2=0;l2<n;l2++) {
      vtmp.y = 0.5*nodes[l2];
      for (l3=0;l3<n;l3++) {
	nodepos[count].x = vtmp.x;
	nodepos[count].y = vtmp.y;
	nodepos[count].z = 0.5*nodes[l3];
	count++;
      }
    }
  }

  // Compute the upward pass matrix
  int n6 = n3*n3;
  double *UpMat = (double *) malloc( n6 *sizeof(double) );
  ComputeWeightsPBC(UpMat, n, alpha, Tkz);

    
  int dof2 = dof.f * dof.s;
  char transa='n';
  double a = 1, b = 0;
  int incr=1;

  // Set the initial scaling factor  
  double scale = 1.0;
        
  vec3 sourcepos1, sourcepos2, origin;
  origin.x = 0, origin.y = 0, origin.z = 0;

  // Point to point interaction
  double Kij1[dof2], Kij2[dof2], Kij3[dof2], Kij4[dof2];
	
  // Ulevel is U^i at the ith level,
  double *Ulevel   = (double *) calloc( n6,  sizeof(double) );
  double *MMresult = (double *) malloc( n6 * sizeof(double) );
        
  // Initialize Ulevel to be identity matrix
  for (k=0; k<n3; k++) 
    Ulevel[k*(n3+1)] = 1.;

  // Allocate and initialize KPBC
  double *KPBC = (double *) calloc( dof2n6, sizeof(double) );

  // double KPBC_j[dof2][n6];
  double *KPBC_j = (double *) malloc( dof2n6 * sizeof(double) );        
          
  // Compute the field values due to periodic sources in different shells
  int ipbc, idof, isource, ifield;
  for (ipbc=1;ipbc<lpbc;ipbc++) {       

    // Compute the jth column that M(j) = KPBC_j * U^i(j)
    for (j=0;j<n3;j++) {
            
      // Initialize KPBC_j
      for (i=0; i<dof2n6; i++)
	KPBC_j[i] = 0.;
            
      // Compute KPBC_j
      for (l=0; l<n3; l++) {
              
	for (l1=-4;l1<5;l1++) {
	  cshift.x = (double)l1;
	  for (l2=-4;l2<5;l2++) {
	    cshift.y = (double)l2;
	    for (l3=-4;l3<5;l3++) {
	      cshift.z = (double)l3;
	      if (abs(l1) > 1 || abs(l2) > 1 || abs(l3) > 1) {
                    
		sourcepos1.x = (cshift.x+nodepos[l].x) * L * scale; 
		sourcepos1.y = (cshift.y+nodepos[l].y) * L * scale; 
		sourcepos1.z = (cshift.z+nodepos[l].z) * L * scale; 
                                            
		sourcepos2.x = sourcepos1.x - nodepos[j].x;
		sourcepos2.y = sourcepos1.y - nodepos[j].y;
		sourcepos2.z = sourcepos1.z - nodepos[j].z;
                      
		kfun(&origin, &sourcepos1, Kij1);
		kfun(&origin, &sourcepos2, Kij2);

		for (k=0; k<n3; k++) {
                                                                   
		  kfun(&nodepos[k], &sourcepos1, Kij3);			
		  kfun(&nodepos[k], &sourcepos2, Kij4);
			
		  for (idof=0; idof < dof2; idof++)
		    KPBC_j[idof*n6 + l*n3 + k] += Kij3[idof] +
		      Kij2[idof] - Kij4[idof] - Kij1[idof];

		}
	      }
	    }
	  }
	}
      }

            
      // Compute M_j
      for (isource=0; isource < dof.s; isource++)
	for (ifield=0; ifield < dof.f; ifield++) {
	  
	  dgemv_(&transa, &n3, &n3, &a, KPBC_j+(isource*dof.f + ifield)*n6, &n3, Ulevel+j*n3, &incr, &a,
		 KPBC+j*dofn3_field*dof.s + ifield + dofn3_field*isource,
		 &dof.f);
	}
            
    }  // j          

          
    // Update 'Ulevel'
    dgemm_(&transa, &transa, &n3, &n3, &n3, &a, Ulevel, &n3,
	   UpMat, &n3, &b, MMresult, &n3);

    Ulevel = Assign(Ulevel, MMresult, n6);
          
    // Increase scaling factor
    scale *= 3.0;		
  } // ipbc
        
  free(UpMat), UpMat=NULL;
  free(Ulevel), Ulevel=NULL;
  free(KPBC_j), KPBC_j=NULL;
  free(MMresult), MMresult=NULL;
        
  //char KPBCfile[50];
  //sprintf(KPBCfile, "Kn%dpbc%d.out", n, lpbc);
  FILE *p_KPBC = fopen(PBCfilename, "w");

  // Write the elastic constants
  for (i=0; i<3; i++)
    fprintf(p_KPBC, "%.10e\n", AnisoParameters->c3[i]);
        
  // Write pbc kernel
  for (i=0; i<dof2n6; i++)
    fprintf(p_KPBC, "%.10e\n", KPBC[i]);
  fclose(p_KPBC);

  free(KPBC);
  printf("Creating pbc file finished.");
}  


/*
 * Function: ComputeWeightsPBC
 * ------------------------------------------------------------------
 * Computes the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) (for PBC calculation each cell has 27 children
 * instead of 8).
 */
static void ComputeSn(vec3 *point, int N, vec3 *Sn, int n, double *Tkz, int grid_type);
void ComputeWeightsPBC(double *UpMat, int n, double alpha, double *Tkz) {

  int i, j, k, l, m, l1, l2, l3, count1, count2;
  vec3 vtmp;

  int n3 = n*n*n;  // n3 = n^3
  int n6 = n3*n3;  // n6 = n^6
  int Nc = 27*n3;  // Number of child Chebyshev nodes
  double nodes[n]; //, vec[n]; 
  vec3 fieldt[Nc]; // Chebyshev-transformed coordinates
  vec3 Sn[n*Nc];
	
  double pi    = M_PI;
  double third = 1.0/3.0; // the parent box is divided into three children in each dimension
	
  // Compute the n Chebyshev nodes of T_n(x)
  for (m=0;m<n;m++)
    nodes[m] = cos(pi*((double)m+0.5)/(double)n);

  // Map all Chebyshev nodes from the children cells to the parent
  k = 0;
  for (i=0;i<27;i++) {
    // Determine the mapping function for the specific child cell
    if (i<9) {
      vtmp.x = -2;			
      if (i<3) 
	vtmp.y = -2;
      else if (i<6)
	vtmp.y = 0;
      else
	vtmp.y = 2;
    } else if (i<18) {
      vtmp.x = 0;
      if (i<12)
	vtmp.y = -2;
      else if (i<15)
	vtmp.y = 0;
      else
	vtmp.y = 2;
    } else {
      vtmp.x = 2;
      if (i<21)
	vtmp.y = -2;
      else if (i<24)
	vtmp.y = 0;
      else
	vtmp.y = 2;
    }
    if (i%3 == 0)
      vtmp.z = -2;
    else if (i%3 == 1)
      vtmp.z = 0;
    else
      vtmp.z = 2;


    // Map all Chebyshev nodes from the children cells to the parent
    // Note: 'nodes' below corresponds to alpha=0
    // scaled child grids: nodes[i]*(1+alpha) +- 2
    // scaled parent box: (1+alpha)*6
    // ihalfL: 1/((1+alpha)*3)
    // so the relative position is: (nodes[i] +- 2/(1+alpha)) / 3
    double L = AdjustBoxSize(1.0, alpha);
    for (l1=0;l1<n;l1++) {
      for (l2=0;l2<n;l2++) {
	for (l3=0;l3<n;l3++) {
	  fieldt[k].x = (nodes[l1] + vtmp.x / L) * third;
	  fieldt[k].y = (nodes[l2] + vtmp.y / L) * third;
	  fieldt[k].z = (nodes[l3] + vtmp.z / L) * third;
	  k++;
	}
      }
    }
  }
    
  // Compute Sc, the mapping function for the field points
  ComputeSn(fieldt, Nc, Sn, n, Tkz, 1); // 1 for Chebyshev
	
  // Compute Sxyz, the weights for the sources
  double *Sxyz = (double *)malloc(n3 * Nc * sizeof(double));
  for (k=0;k<Nc;k++) {
    l = 0;
    for (l1=0;l1<n;l1++) {
      for (l2=0;l2<n;l2++) {
	for (l3=0;l3<n;l3++) {
	  Sxyz[l*Nc+k] = Sn[l1*Nc+k].x*Sn[l2*Nc+k].y*Sn[l3*Nc+k].z;
	  l++;
	}
      }
    }
  }


  // Accumulate the weights into a single weight matrix W
  double *W = (double *)calloc(n6, sizeof(double));
	
  for (i=0;i<27;i++) {
    count1 = 0;
    count2 = i*n3;
    for (j=0;j<n3;j++) {
      for (k=0;k<n3;k++) {
	W[count1] += Sxyz[k*Nc+count2]; // Notice the transpose here
	count1++;
      }
      count2++;
    }
  }
	
  UpMat = Assign(UpMat, W, n6);
   
  free(Sxyz);
  free(W);	
}


// summation using Clenshaw's recurrence relation
static double ClenshawSum(int n, double x, double *Tk) {
  int j;
  double d0, d1, d2;
  d0 = d1 = 0;
  for (j = n - 1; j > 0; j--) {
    d2 = d0;
    d0 = 2.0 * x * d0 - d1 + Tk[j];
    d1 = d2;
  }
  return x * d0 - d1 + 0.5 * Tk[0];
}

void ComputeSn(vec3 *point, int N, vec3 *Sn, int n, double *Tkz, int grid_type) {
     
  int i, k, m;

	    
    double pfac = 2./n;	    
    double *Tkz_m;

    // Loop over Chebyshev node m
    for (m=0;m<n;m++) {
      k    = m*N;
      Tkz_m = Tkz + m*n;
      for (i=0;i<N;i++) {
	// Compute S_n for each direction using Clenshaw
	Sn[k+i].x = pfac * ClenshawSum(n, point[i].x, Tkz_m);
	Sn[k+i].y = pfac * ClenshawSum(n, point[i].y, Tkz_m);
	Sn[k+i].z = pfac * ClenshawSum(n, point[i].z, Tkz_m);
      }
    }

} // end function

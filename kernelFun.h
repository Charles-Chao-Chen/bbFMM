/*
 * This file provides different kernel functions to be used in the fast
 * multipole method
 */

#ifndef _kernelfun_h
#define _kernelfun_h


#include "AnisoFunctions.h" // real8, real8_complex
#include "utility.h"


// define function pointers
typedef void (*kfun_t) (vec3*, vec3*, double*);

typedef struct {
  char name[50];
  double homogen;
  int symm;
  //int2 dof;
  kfun_t kfun;
} kernel_t;

// example kernels are:
// kernel_t Lapforce   = {"lapforce",  2.0, -1, &LapforceFun};
// kernel_t Laplacian  = {"laplacian", 1.0,  1, &LaplacianFun};
// kernel_t Gaussian   = {"gaussian",  0.0,  1, &GaussianFun};
// kernel_t Polynomial = {"poly",     -3.0,  0, &PolyFun};
// kernel_t OneOverR4  = {"1overr4",   4.0,  1, &OneOverR4Fun};
// kernel_t DDANISO    = {"Aniso",     2.0, -1, &Fun};


void PolyFun (vec3* fieldpos, vec3* sourcepos, double *K);

void GaussianFun (vec3* fieldpos, vec3* sourcepos, double *K);

void LaplacianFun (vec3* fieldpos, vec3* sourcepos, double *K);
  
// test PBC implementation 
void LapforceFun (vec3* fieldpos, vec3* sourcepos, double *K);


// initialize anisotropic parameters
void AnisoInit();


// Anisotropic kernel function in dislocation dynamics
void AnisoFun (vec3* fieldpos, vec3* sourcepos, double *K);


// clean up anisotropic parameters
void AnisoClean();



// Add interaction from PBC using FMM
void DoCorrectionTable(double *ChebyshevWeightSource, double *ChebyshevWeightField, int n, int2 dof, double Len, double alpha, int lpbc, kernel_t kernel, double *Tkz);


void CreateTableCorrection(kfun_t kfun, int n, int2 dof, double alpha, int lpbc, double *Tkz, char* PBCfilename);

  
void ComputeWeightsPBC(double *UpMat, int n, double alpha, double *Tkz);

/*
double* MeanStressFmm(const double *ChebyshevWeightField, const int n, const int dof, const double *Tkz);
*/
  
// Discretize segments into points using Gauss quadrature
void Segment2Point(segT *segment, int segmentSize, double *burg, int
		   numGauss, vec3 *source, double *sourceCharge);


//#ifdef TENSOR
/*
  void DoCorrectionTable(vec3 *field, int Nf, vec3 *source, 
  int Ns, double *sourceCharge, int n, deg2
  *dof, double Len, double elasConst[3][3][3][3], 
  int lpbc, double *fieldStress);
*/
/*
  #elif LINEINT
  void DoTableCorrection(vec3 *field, int Nf, point2 *segment, int
  nSegment, int numGauss, double *burg, 
  int n, deg2 *dof, double Len, double elasConst[3][3][3][3], 
  int lpbc, double *fieldStress);

  #endif
*/
/*
  void FMMPBC(vec3 *field, int Nf, vec3 *source, int Ns, double *sourceCharge, 
  int n, deg2 *dof, double Len, double *Tkz, double elasConst[3][3][3][3], 
  int lpbc, double *fieldStress);
*/

/*
 * Function: ComputePeriodicKernel
 * ---------------------------------------------------------------------
 * Forms the matrix that describes the interactions of the computational
 * cell with its periodic images up to lpbc shells.
 */
//void ComputePeriodicKernel(double *KPBC, double *UpMat, 
//      double elasConst[3][3][3][3], 
//    double L, int n, deg2 *dof,int lpbc);
//void CreateTableCorrection(kfun_t kfun, int n, deg2 *dof,int lpbc, double *Tkz, int grid_type);


// parameters for evaluating anisotropic kernel
typedef struct _paraAniso{

  int sizeReFq, sizeImFq, qMax;
  int *ind1Re, *ind2Re, *ind1Im, *ind2Im;

  double elasA;
  real8_complex *e12;
  real8 (*CdGdx)[6];
  double *c3, *e3, (*C633)[3][3]; //, *C;
  double *ReFqout, *ImFqout;

  char *f_Ind1Re, *f_Ind2Re, *f_Ind1Im, *f_Ind2Im, *f_ReFqout,
    *f_ImFqout;    
    
} paraAniso;


#endif



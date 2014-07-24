#ifndef _BBFMM_H
#define _BBFMM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>  // false
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <fftw3.h>

#include "AnisoFunctions.h"
#include "kernelFun.h"
#include "utility.h"

#ifdef __cplusplus
extern "C" {
#endif
  

  
  typedef enum {
    UNIFORM,
    CHEBYSHEV,
  } gridT;
  

  typedef enum {
    POINT,
    SEG,
  } srcT;

  
  typedef struct {
    srcT src_t;
    segT *seg; // segment source
    vec3 *source; // point source
    double *sigma;  // source intensity
    vec3 *target; // target points
    int n;  // interpolation order
    int Ns; // source size
    int Nf; // target size
    double boxL;  // box length
    double boxA;  // box adjustment parameter for segment source
    int fmm_lvl;  // fmm tree level
    int pbc_lvl;  // fmm pbc level
    gridT grid_t; // interpolation grid type
  } fmm_input;

  
  /*
   * Function: Timer
   * ----------------------------------------
   * Returns the time in seconds.
   *
   */
  typedef double timeType;
  timeType Timer(void);
#define evaltime(timeval_time) (double)timeval_time.tv_sec	\
  + (double)timeval_time.tv_usec*1e-6

  /* Struct: nodeT
   * -------------------------------------------------------------------
   * This struct is a node in the octree.
   */
  typedef struct _nodeT {
    struct _nodeT *leaves[8];
    struct _nodeT *parent;
    struct _nodeT *neighbors[27];
    struct _nodeT *interaction[189];
    vec3 center, cshiftneigh[27], cshiftinter[189];
    double length;
    double *fieldval, *sourceval, *sourcefre, *proxysval; 
    int *fieldlist, *sourcelist;
    int Nf, Ns, ineigh, iinter;
  } nodeT;


  /*
    typedef struct _dof_struct {
    int s; // Size of source vector
    int f; // Size of field vector
    } doft; // Type of variable dof
  */

  /* Struct: fmmparam
   * -------------------------------------------------------------------
   * This struct stores all of the FMM parameters.
   */
  
  typedef struct _fmmparam {
    int Ns;       // Number of sources
    int Nf;       // Number of field points
    dof2 dof;        // Number of degrees of freedom
    double L;        // Length of one side of simulation cube
    int n;        // Number of Chebyshev nodes in one direction
    int levels;        // Maximum number of levels in octree
    int PBClevels;        // Number of PBC levels (images = 27^PBClevels)
    int PBCshells;        // Number of PBC shells (for direct calculation)
    int precomp;        // Turn on (1) or off (0) pre-computation
    dof2 cutoff;       // Number of singular values to keep
    double homogen;        // Order of homogeneity of kernel
    char filesval[80];     // File name for storing singular values
    char filesvec[80];     // File name for storing singular vectors
  } fmmparam;



  /*
   * Function: bbfmm
   * -----------------------------------------------------------------
   * Given the source and field point locations, strength of the sources,
   * number of field points and sources, length of the computation cell, and 
   * the number of Chebyshev nodes, the field is computed.
   */

#ifdef LINEINT
  void bbfmm(vec3 *field, segT *segment, double *burg, int Nf,
	     int Ns, dof2 dof, double boxLen, double alpha, int
	     level, int n, kernel_t kernel, double epsilon,
	     int numgau, double *phi, int grid_type, int lpbc);

#elif TENSOR
double* bbfmm( vec3 *field, int Nf, vec3 *source, int Ns, double *q,
	       dof2 dof,
	       double box_len, double alpha, int n, gridT grid,
	       kernel_t kernel, int level_tree, int level_pbc,
	       double epsilon );

#endif


  /*
   * Function: FMMSetup
   * -----------------------------------------------------------------
   * Prepare for the FMM calculation by setting the parameters, computing
   * the weight matrices, pre-computing the SVD (if necessary), reading 
   * in the necessary matrices, and building the FMM hierarchy.
   */
void FMMSetup(nodeT **A, double *Tkz, int *Ktable, double
	      *Kweights, double *Cweights, double boxLen, double
	      alpha, dof2 *cutoff, int n,
	      kernel_t kernel, double epsilon, dof2
	      dof, int Ns, int treeLevel, char *Kmat, char *Umat,
	      char *Vmat, int grid_type);
  

bool precompute_files_usable( char *Kmat, char *Umat, char *Vmat,
			      double homogen, double boxLen,
			      int treeLevel, int grid_type );


  void generate_precompute_files(double boxLen, int treeLevel, int n, dof2 dof, kernel_t kernel, char *Kmat, char *Umat, char *Vmat, double alpha, double *Kweights, double epsilon, int grid_type);

     
  void compute_m2l_operator (int n, dof2 dof, kernel_t kernel, char *Kmat, char *Umat, char *Vmat, double l, double alpha, double *Kweights, double epsilon, int grid_type);
     
  /*
   * Function: FMMReadMatrices
   * ------------------------------------------------------------------
   * Read in the kernel interaction matrix M and the matrix of singular
   * vectors U.
   */
  void FMMReadMatrices(double **K, double **U, double **VT, dof2 cutoff,
		       int n, dof2 dof, char *Kmat, char *Umat,
		       char *Vmat, int treeLevel, double homogen,
		       int grid_type);



  /*
   * Function: FMMDistribute
   * ------------------------------------------------------------------
   * Distribute the field points and sources to the appropriate location
   * in the FMM hiearchy and sets up the interaction list.
   */
  void FMMDistribute(nodeT **A, vec3 *field, vec3 *source,
		     int Nf, int Ns, int level);

  /*
   * Function: FMMCompute
   * ------------------------------------------------------------------
   * Computes the field using BBFMM.
   */

#ifdef LINEINT
  void FMMCompute(nodeT **A, vec3 *field, int Nf, segT *segment, 
		  vec3 * midpoint, double *burg, double *K, double *U, 
		  double *VT, double *Tkz, int *Ktable, double *Kweights,
		  double *Cweights, dof2 cutoff, int n, 
		  dof2 dof, int numgau, double alpha, double *phi,
		  int grid_type, double boxLen, int lpbc,
		  kernel_t kernel);

#elif TENSOR
  void FMMCompute(nodeT **A, vec3 *field, int Nf, vec3 *source, 
		  double *q, double *K, double *U, 
		  double *VT, double *Tkz, int *Ktable, double *Kweights,
		  double *Cweights, dof2 cutoff, int n, dof2 dof,
		  double alpha, double *phi, int
		  grid_type, double boxLen, int lpbc,
		  kernel_t kernel);
#endif

  /*
   * Function: FMMCleanup
   * ------------------------------------------------------------------
   * Cleans up after FMM calculation.
   */
  void FMMCleanup(nodeT *A);

  /*
   * Function: SetParam
   * ------------------------------------------------------------------
   * Read in the user-specified options from the options file.
   * 
   */
  int SetParam(char *options, fmmparam *param);


  /* This function computes the Gaussian points along the segments in 
   * a leaf box with their quadrature weights
   */
  void LineIntegral(segT *segment, vec3 *midpoint, int Ns,
		    int *sourcelist, int num, double *gpoint, double
		    L, vec3 scenter, vec3 *sourcet, double *xi);


  /*
   * Function: EvaluateKernel
   * -------------------------------------------------------------------
   * Evaluates the kernel given a source and a field point.
   */
  void EvaluateKernel(vec3 fieldpos, vec3 sourcepos, paraAniso
		      *AnisoParameters,
		      double *K, dof2 dof);

  /*
   * Function: EvaluateKernelCell
   * -------------------------------------------------------------------
   * Evaluates the kernel for interactions between a pair of cells.
   */
  void EvaluateKernelCell(vec3 *fieldpos, vec3 *sourcepos, 
			  int Nf, int Ns, dof2 dof, kfun_t kfun, double *kernel);

  /*
   * Function: EvaluateField
   * -------------------------------------------------------------------
   * Evaluates the kernel for interactions between a pair of cells.
   */
  void EvaluateField(vec3 *fieldpos, vec3 *sourcepos, 
		     double *q, int Nf, int Ns, dof2 dof, 
		     kfun_t kfun, double *fieldval);

  /*
   * Function: ComputeWeights
   * ------------------------------------------------------------------
   * Computes the weights for the Chebyshev nodes of all children cells
   * (identical for all cells and all levels so just compute once and
   * store in memory) and set up the lookup table.
   */
  void ComputeWeights(double *Tkz, int *Ktable, double *Kweights, 
		      double *Cweights, int n, double alpha, int
		      grid_type);

  
  //void ComputeWeightsPBC(double *UpMat, int n, double *Tkz, int grid_type);


  /*
   * Function: ComputeKernelSVD
   * -------------------------------------------------------------------
   * Computes the kernel for 316n^6 interactions between Chebyshev nodes
   * and then computes the SVD of the kernel matrix.
   */
  void ComputeKernelSVD(double *Kweights, int n, kernel_t kernel,
			double epsilon, dof2 dof,
			char *Kmat, char *Umat, char *Vmat, double
			alpha, double boxLen);


  /*
   * Function: ComputeKernelUniformGrid
   * ------------------------------------------------------------------
   * Computes the kernel for 316(2n-1)^3 interactions between Uniform
   * Grid nodes. Does not compute SVD.
   */
  void ComputeKernelUniformGrid(int n,dof2 dof, kfun_t kfun, char
				*Kmat, double alpha, double len);
 

  /*
   * Function: NewNode
   * ------------------------------------------------------------------
   * Dynamically allocates space for a new node A in the octree and
   * initializes the quantities in the node.
   *
   */
  void NewNode(nodeT **A, vec3 center, double L, int n);

  /*
   * Function: BuildFMMHierarchy
   * ------------------------------------------------------------------
   * Builds the FMM hierarchy with l levels.
   *
   */
  void BuildFMMHierarchy(nodeT **A, int l, int n);

  /*
   * Function: DistributeFieldPoints
   * -------------------------------------------------------------------
   * Distributes all of the field points to the appropriate cell at the 
   * bottom of the octree.
   */
  void DistributeFieldPoints(nodeT **A, vec3 *field, int *fieldlist, int levels);

  /*
   * Function: DistributeSources
   * -------------------------------------------------------------------
   * Distributes all of the sources to the appropriate cell at the 
   * bottom of the octree.
   */
  void DistributeSources(nodeT **A, vec3 *segment, 
			 int *sourcelist, int levels);

  /*
   * Function: InteractionList
   * -----------------------------------------------------------------
   * For each node at each level of the octree, the interaction list
   * is determined.
   */
  void InteractionList(nodeT **A, int levels);


  double AdjustBoxSize(double Len, double alpha);

  /*
   * Function: UpwardPass
   * -------------------------------------------------------------------
   * Gathers the sources from the children cells and determines the strength
   * of pseudo-charges located at the Chebyshev nodes of the parent cell.
   * (upward pass of BBFMM)
   */

#ifdef LINEINT
  void UpwardPass(nodeT **A, segT *segment, vec3 *midpoint,
		  double *Tkz, double *burg, double *Kweights, 
		  int n, int dof, int numgau, double *gpoint,
		  double *gweight, double alpha, int grid_type);
		  //, double boxLen, double homogen, int lpbc, kfun_t kfun);

#elif TENSOR
  void UpwardPass(nodeT **A, vec3 *source, double *q,
		  double *Cweights, double *Tkz,
		  dof2 cutoff, int n, int dof, int grid_type);
#endif

  void Moment2Moment(int n, double *r, double *Schild, double *SS, 
		     int dof, double *Cweights);

  void Particle2Moment(nodeT **A, segT *segment, double *burg, 
		       vec3 *midpoint, int n, double *Tkz, int dof, 
		       int numgau, double *gpoint, double *gweight,
		       double alpha, int grid_type);

  /*
   * Function: FMMInteraction
   * ---------------------------------------------------------------------
   * At each level of the octree the interaction between well-separated cells of 
   * field and source Chebyshev nodes is computed using a SVD lookup table.
   */
  void FMMInteraction(nodeT **A, double *E, int *Ktable, double *U, 
		      double *VT, double *Kweights, int n, dof2 dof,
		      dof2 cutoff, double homogen, int treeLevel, int
		      grid_type);


  void Moment2Local(int n, double *R, double *cell_mpCoeff, 
		    double *FFCoeff, double *E, int *Ktable, dof2
		    dof, dof2 cutoff, double *VT, double *Kweights,
		    int grid_type);


  /*
   * Function: FMMInteractionPBC
   * ---------------------------------------------------------------------
   * Add the interactions due to the periodic images.
   */
  //void FMMInteractionPBC(nodeT **A, double *KPBC, int n, dof2 dof);

  /*
   * Function: DownwardPass
   * -------------------------------------------------------------------
   * Distributes the field from the parent cell to the children cells 
   * using Chebyshev interpolation. (downward pass of BBFMM)
   */

#ifdef LINEINT
  void DownwardPass(nodeT **A, vec3 *field, double *Cweights, 
		    double *Tkz, int n, int dof, double alpha,
		    double *phi, int grid_type);

#elif TENSOR
  void DownwardPass(nodeT **A, vec3 *field, vec3 *source, 
		    double *Cweights, double *Tkz, double *q,
		    dof2 cutoff, int n, dof2 dof,
		    double alpha, kfun_t kfun,
		    double *phi, int grid_type);
                                                 
#endif


  void Local2Local(int n, double *r, double *F, double *Fchild, 
		   int dof, double *Cweights, int grid_type);

  void Local2Particle(nodeT **A, vec3 *field, double *Tkz, int n,
		      int dof, double alpha, double *phi, int
		      grid_type);

  /*
   * Function: FreeNode
   * ------------------------------------------------------------------
   * Frees up space associated with node A.
   *
   */
  void FreeNode(nodeT *A);


  void GetPosition(int n, int idx, double *fieldpos, double *sourcepos,
		   double *nodepos);

  void FrequencyProduct(int N, double *Afre, double *xfre, double *res);



void AddPBC(int lpbc, double *M, double **L, double *phi, int Nf, int n, dof2 dofst, double boxLen, double alpha, double *Tkz, kernel_t kernel, int grid_type );

void GridPos1d(double alpha, double len, int n, int grid_type, double *nodes);
void GridPos3d(const vec3 center, const double alpha, const double L, const int n, const gridT grid, vec3 *P);

void Anterpolate(double *In, int in_gridtype, double *Out, int out_gridtype, int n, int dof, double *Tkz);
void Interpolate(double *In, int in_gridtype, double *Out, int out_gridtype, int n, int dof, double *Tkz);

void Create_lookup_table(int *Ktable);
void ComputeSn(vec3 *point, int N, vec3 *Sn, int n, double *Tkz, int grid_type);
void ComputeTkz(double *T, int n);

inline bool Is_well_separated(vec3 p, double L);

  /*
inline void InBoxCheck( double x ) {
  if (fabs(x) > 1)
    printf("Outside Box!\n");
}
*/

  double* MeanStressFmm(const double *ChebyshevWeightField, const int n, const int dof, const double *Tkz);

  
#ifdef __cplusplus
}
#endif
  
#endif

#ifndef _vec3_h
#define _vec3_h
  
#include <stdbool.h>

// 3d point
typedef struct {
  double x,y,z;
} vec3;


// 3d segment
typedef struct {
  vec3 p1;
  vec3 p2;
} segT;


/* Struct: deg2
 * -------------------------------------------------------------------
 * The FMM can handle cases where the input are output variables
 * are vectors (not just scalars).
 * This struct stores information about the size of the input
 * vector (e.g., 1, 3, 9, etc) and the output vector (e.g., 1, 3, 6).
 * This is different from the number of source and field points.
 */
// source and field degree
typedef struct {
  int s; // Size of source vector
  int f; // Size of target vector
} dof2;


/* Uniform random number generator */
#define frand(xmin,xmax) ((double)xmin+(double)(xmax-xmin)*rand()/	\
			  (double)RAND_MAX) 


// x+=y
inline double* Add(double *x, const double *y, const int n);

// x-=y
inline double* Subtract(double *x, const int n1, const double *y, const int n2);

// x*=a
inline double* Multiply(double *x, const double a, const int n);
  
// x<-y
inline double* Assign(double *x, const double *y, const int n);

// set all elements to zero
void Zero(double *x, const int N);

// L2 norm
inline double Norm(const double *x, const int n);


inline bool equal_double(double, double);

inline bool equal_vec3(vec3*, vec3*);


void print_array( const double *arr, const int N, const char *arr_name );
	  

// create the directory if not exist
inline void create_directory( char *name);


#ifdef __cplusplus
extern "C" {
#endif

  // Declaration for BLAS dot product
  double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);

  // Declaration for BLAS matrix-matrix multiply
  void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	      double *alpha, double *A, int *lda, double *B,
	      int *ldb, double *beta, double *C, int *ldc);

  // Declaration for BLAS matrix-vector multiply
  void dgemv_(char *trans, int *m, int *n, double *alpha, double *A,
	      int *lda, double *x, int *incx, double *beta, 
	      double *y, int *incy);

  // Declaration for LAPACK's SVD routine
  void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *A,
	       int *lda, double *S, double *U, int *ldu, double *VT,
	       int *ldvt, double *work, int *lwork, int *info);
  
#ifdef __cplusplus
}
#endif


/*
 * Function: ComputeError
 * ------------------------------------------------------------------
 * Computes the L2 relative error. Note: x is changed here.
 */
double ComputeL2Err(double *x, double *y, int N);


/*
 * Function: SetSources  -> to be replaced  x' = sources, field = x
 * -------------------------------------------------------------------
 * Distributes an equal number of positive and negative charges uniformly
 * in the simulation cell ([-0.5*L,0.5*L])^3 and take these same locations
 * as the field points.
 */
void SetSources(vec3 *field, int Nf, vec3 *source, int Ns, double *q,
		int dof, double L);

void SetSourcesCase1(vec3 *field, vec3 *source, double *q,
		     int N, dof2 * dof, double L, double r);

//void SetSourcesCase2(vec3 *field, vec3 *source, double *q,
//		       int N, dof2 * dof, double L, vec3 *scent, double *rrr);

void SetSourcesCase3(vec3 *field, vec3 *source, double *q,
		     int N, dof2 * dof, double L, vec3 *scent, double *rrr);


#endif


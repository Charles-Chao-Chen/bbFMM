#ifndef _utility_h
#define _utility_h
  
#include <stdbool.h>


#define debugging_enabled 1

#define DEBUG(x) do {				\
    if (debugging_enabled) {			\
      printf("[%s: %i] %s\n",			\
	     __FILE__, __LINE__, x);		\
    }						\
  } while (0)					\
    

#define READ_CHECK( callReturn, num ) do {		\
    if (callReturn != num) {				\
      printf("Read error in file '%s' at line %i.\n"	\
	     , __FILE__, __LINE__);			\
      exit(1);						\
    }							\
  } while(0)						\


#define FILE_CHECK( fptr, fname ) do {				\
    if (fptr == NULL) {						\
      printf("Can't open file '%s' in '%s' at line %i.\n"	\
	     , fname, __FILE__,  __LINE__);			\
      assert(false);						\
    }								\
  } while(0)							\



// global varibale for Gauss quadrature
extern double *GAUSSP;
extern double *GAUSSW;


/* Function: Timer
 * ----------------------------------------
 * Returns the time in seconds.
 */
typedef double timeType;
timeType Timer(void);
#define evaltime(timeval_time) (double)timeval_time.tv_sec	\
  + (double)timeval_time.tv_usec*1e-6


// 3d point
typedef struct {
  double x,y,z;
} vec3;


// 3d segment
typedef struct {
  vec3 p_beg;
  vec3 p_end;
} segT;


/* Struct: int2
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
} int2;


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

// initialize Gauss quadrature points and weights
void InitGaussQuadrature(int N);
void CleanGaussQuadrature();

#ifdef __cplusplus
extern "C" {
#endif

  /* Note that blas/lapack routines are Fortran-style:
   *   (1) Pass variables by address, not by value.
   *   (2) Store your data in Fortran style, that is, column-major
   *   rather than row-major order.
   */
  
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
		     int N, int2 * dof, double L, double r);

//void SetSourcesCase2(vec3 *field, vec3 *source, double *q,
//		       int N, int2 * dof, double L, vec3 *scent, double *rrr);

void SetSourcesCase3(vec3 *field, vec3 *source, double *q,
		     int N, int2 * dof, double L, vec3 *scent, double *rrr);


#endif


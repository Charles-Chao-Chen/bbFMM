#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>       // for 'struct timeval', 'gettimeofday()' and 'evaltime()'
#include <sys/types.h>      // create directory
#include <sys/stat.h>       // create directory
#include <math.h>
#include <time.h>           // for srand()

#include "AnisoFunctions.h" // for 'gqwp()'
#include "utility.h"

#define EQUAL_DOUBLE_EPS 1e-6


/* Function: Timer
 * ----------------------------------------
 * Returns the time in seconds.
 */
extern timeType Timer(void) {
  struct timeval timeval_time;
  gettimeofday(&timeval_time,NULL);
  return evaltime(timeval_time);
}


double* Add(double *x, const double *y, const int n) {

  int i;
  for (i=0; i<n; i++) {
    x[i] += y[i];
  }
  return x;
}


double* Subtract( double *x, const int n1, const double *y, const int n2 ) {

  assert( n1%n2 == 0);
  int N = n1/n2, count = 0, i, j;
  for (i=0; i<N; i++)
    for (j=0; j<n2; j++) {
      x[count] -= y[j];
      count ++;
    }
	    
  return x;
}

double* Multiply(double *x, const double a, const int n) {

  int i;
  for (i=0; i<n; i++)
    x[i] *= a;

  return x;
}


double* Assign( double *x, const double *y, const int n ) {
  int i;
  for (i=0; i<n; i++)
    x[i] = y[i];
  return x;
}


void Zero(double *x, const int N) {
  int i;
  for (i=0; i<N; i++)
    x[i] = 0;
}

void create_directory( char *name ) {
    struct stat st = {0};
    if (stat(name, &st) == -1) {
      printf("Create '%s' dirctory.\n", name);
      mkdir(name, 0775);
    }
}


//int     NGAUSS;
double *GAUSSP;
double *GAUSSW;

//void SetGaussOrder(int N) { NGAUSS = N;}

void InitGaussQuadrature(int N) {
  GAUSSP = malloc(N*sizeof(double));
  GAUSSW = malloc(N*sizeof(double));
  gqwp(N, GAUSSP, GAUSSW);
}


void CleanGaussQuadrature() {
  free(GAUSSP);
  free(GAUSSW);
}


void print_array( const double *arr, const int N, const char *arr_name ) {

  printf("%s:\n", arr_name);
  int i;
  for (i=0; i<N; i++)
    printf("%f\t", arr[i]);
  printf("\n");
}


bool equal_double(double a, double b) {
  return fabs(a - b) < EQUAL_DOUBLE_EPS;
}


bool equal_vec3(vec3* p1, vec3* p2) {

  return
    equal_double( p1->x, p2->x ) &&
    equal_double( p1->y, p2->y ) &&
    equal_double( p1->z, p2->z );
}


double Norm(const double *x, const int n) {
  double sum = 0;
  int i;
  for (i=0; i<n; i++)
    sum += x[i]*x[i];
  return sqrt(sum);
}


double ComputeL2Err(double *x, double *y, int N) {

  x = Subtract( x, N, y, N );
  double s1 = Norm(x, N), s2 = Norm(y, N);
  //printf("norm(x): %.2f, norm(y): %.2f\n", s1, s2);
  return s1/s2;
}


/*
 * Function: SetSources  -> to be replaced  x' = sources, field = x
 * -------------------------------------------------------------------
 * Distributes an equal number of positive and negative charges uniformly
 * in the simulation cell ([-0.5*L,0.5*L])^3 and take these same locations
 * as the field points.
 */
void SetSources(vec3 *field, int Nf, vec3 *source, int Ns, double *q,
		int dof, double L) {

  //srand( time(NULL) );
  
  int i, j, k=0;
	
  // Distributes the sources randomly uniformly in a cubic cell
  for (i=0;i<Ns;i++) {
    for (j=0; j<dof; j++)
      q[k++] = frand(0,1);
		
    source[i].x = frand(-0.5,0.5)*L;
    source[i].y = frand(-0.5,0.5)*L;
    source[i].z = frand(-0.5,0.5)*L;
  }
	
  // Randomly set field points
  for (i=0;i<Nf;i++) {
    field[i].x = frand(-0.5,0.5)*L;
    field[i].y = frand(-0.5,0.5)*L;
    field[i].z = frand(-0.5,0.5)*L;
  }
}

// Points distributed outside a cube of radius r and inside the unit cube
void SetSourcesCase1(vec3 *field, vec3 *source, double *q,
		     int N, int2 * dof, double L, double r) {
	
  int i, j, k;
  double tmp;
	
  k = 0;
	
  double rr = 0.5 - r;
  for (i=0;i<N;i++) {
    if (i%2 == 0)
      tmp = 1;
    else
      tmp = -1;
		
    for (j=0; j<dof->s; j++) {
      q[k] = tmp;
      k++;
    }
		
    source[i].x = frand(-rr,rr)*L;
    source[i].x = source[i].x > 0 ? source[i].x+r : source[i].x-r;
    source[i].y = frand(-rr,rr)*L;
    source[i].y = source[i].y > 0 ? source[i].y+r : source[i].y-r;
    source[i].z = frand(-rr,rr)*L;
    source[i].z = source[i].z > 0 ? source[i].z+r : source[i].z-r;
  }
	
  // Takes the source locations as the field points
  for (i=0;i<N;i++) {
    field[i].x = source[i].x;
    field[i].y = source[i].y;
    field[i].z = source[i].z;
  }
}

/*
  void SetSourcesCase2(vec3 *field, vec3 *source, double *q,
  int N, int2 * dof, double L, vec3 *scent, double *rrr) {
 
  int i, j, k;
  double tmp;
 
  // Distributes the sources uniformly in the cubic cell
  k = 0;
  for (i=0;i<N;i++) {
  if (i%2 == 0)
  tmp = 1;
  else
  tmp = -1;
 
  for (j=0;j<dof->s;j++) {
  q[k] = tmp;
  k++;
  }
 
  int ind = i%8;
  source[i].x = frand(scent[ind].x-rrr[ind],scent[ind].x + rrr[ind])*L;
  source[i].y = frand(scent[ind].y-rrr[ind],scent[ind].y + rrr[ind])*L;
  source[i].z = frand(scent[ind].z-rrr[ind],scent[ind].z + rrr[ind])*L;
  }
 
  // Takes the source locations as the field points
  for (i=0;i<N;i++) {
  field[i].x = source[i].x;
  field[i].y = source[i].y;
  field[i].z = source[i].z;
  }
  }
*/

void SetSourcesCase3(vec3 *field, vec3 *source, double *q,
		     int N, int2 * dof, double L, vec3 *scent, double *rrr) {
	
  int i, j, k;
  double tmp;
	
  k = 0;
  for (i=0;i<N;i++) {
    if (i%2 == 0)
      tmp = 1;
    else
      tmp = -1;
		
    for (j=0; j<dof->s; j++) {
      q[k] = tmp;
      k++;
    }
		
    // Points randomly distributed in 8 subcubes in level 1 of the hierarchy tree  
    int ind = i%8;
    source[i].x = frand(scent[ind].x-rrr[ind],scent[ind].x + rrr[ind])*L;
    source[i].y = frand(scent[ind].y-rrr[ind],scent[ind].y + rrr[ind])*L;
    source[i].z = frand(scent[ind].z-rrr[ind],scent[ind].z + rrr[ind])*L;
  }
	
  // Take the source locations as the field points
  for (i=0;i<N;i++) {
    field[i].x = source[i].x;
    field[i].y = source[i].y;
    field[i].z = source[i].z;
  }
}

#include "utility.h"
#include "kernelFun.h"


/*
 * Function DirectCalc3D
 * -------------------------------------------------------------------
 * Computes the potential at the first field point and returns 
 * the result in phi.
 */
double* directCalc( vec3 *field, int Nf, vec3 *source, int Ns, double *q,
		    dof2 dof, double box_len, kfun_t kfun, int level_pbc );

/*
void DirectCalc3D(vec3 *field, int Nf, vec3 *source, double
		  *sourceCharge, int sourceSize, dof2 dof, int levelpbc,
		  double boxLen, kfun_t kernel, double
		  *stressField);
*/

// direct result without lpbc=1
void DirWo1 (double *phi, int N, int lpbc);

// direct result of pbc
void DirPBC (double *phi, int N, int lpbc);

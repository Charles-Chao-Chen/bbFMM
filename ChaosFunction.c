#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "AnisoFunctions.h"

void ChaoStress(double p1x, double p1y, double p1z,
		double p2x, double p2y, double p2z,      
		int qMax,
		int *indexList1Real, int *indexList2Real,
		int *indexList1Imag, int *indexList2Imag,
		real8 *FqReal, real8 *FqImag,
		int FqRealNumElem, int FqImagNumElem,
		real8_complex e12[3],real8 e3[3], real8 C633[6][3][3],
		double TrapdGdx[3][6],double K[54])
{
  /*
    For the anisotropic kernel, K has size dof.f x dof.s (6 x 9).
    K is stored in column-major format.
  */
	
  //printf("e12=(%f, %f, %f), ", e12[0], e12[1], e12[2]);
  //printf("e3=(%f, %f, %f), ", e3[0], e3[1], e3[2]);
  //assert(e3[2] == 1);
  
  /* The integral version (numerical quadrature) for computing dG/dx  */
  DerivativeOfGreensFunction3(p1x, p1y, p1z, 
			      p2x, p2y, p2z, 
			      qMax, 			  
			      indexList1Real,indexList2Real,
			      indexList1Imag,indexList2Imag,
			      FqReal,FqImag,
			      FqRealNumElem,FqImagNumElem,
			      e12,e3,TrapdGdx);
  int v,w,n,js;
  int alpha;
  double D[3][3][3];

  for (v=0; v<3; v++) 
    for (w=0; w<3; w++) 		  
      for (n=0; n<3; n++) { 
	D[v][w][n] = 0.0;
	for (alpha=0;alpha<6;alpha++) {	  
	  D[v][w][n] += C633[alpha][w][n]*TrapdGdx[v][alpha];
	}
      }


  // (j,s) and (w,r) indices
  for (w=0; w<3; ++w) 			
    for (js=0; js<6; ++js) {
      K[js+w*6] = 0.0;	      
      K[js+w*6+18] = 0.0;	      
      K[js+w*6+2*18] = 0.0;	      
      for (v=0; v<3; v++) {
	K[js+w*6]      += C633[js][v][2] * D[v][w][1]-C633[js][v][1] * D[v][w][2];
	K[js+w*6+18]   += C633[js][v][0] * D[v][w][2]-C633[js][v][2] * D[v][w][0];
	K[js+w*6+2*18] += C633[js][v][1] * D[v][w][0]-C633[js][v][0] * D[v][w][1];
      }
    }
}

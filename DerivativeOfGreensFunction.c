#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <time.h>

#include "AnisoFunctions.h"

/*---------------------------------------------------------------------------
 *
 *      Function:    DerivativeOfGreensFunction
 *      Description: Derivative of the Green's function using numerical 
 *                   integration formula in [0, pi] evaluated at R = x2-x1.
 *
 *      Parameters:
 *          IN:     nMax
 *          IN:     x1   : A field point
 *          IN:     x2   : A point on the segment
 *          IN:     C    : the elastic tensor in the form of a 3x3x3x3 matrix
 *          OUT:    dGdx : the derivative of the Green's function. The last index
 *                         is the index for the derivative dG_ij dx_k = dGdx[3][3][3]
 *
 *-------------------------------------------------------------------------*/
void DerivativeOfGreensFunction(int nMax, real8 C[3][3][3][3],
				real8 p1x, real8 p1y, real8 p1z,
				real8 p2x, real8 p2y, real8 p2z,
				real8 dGdx2[3][6])
{
        int   i, j, k, l, m, n, p;
        int   vecMinIndex;
        real8 tmp1, vec[3], dX, dY, dZ, temp, invL;
        real8 xHat[3], yHat[3];
	real8 dGdx[3][3][3];


        dX = p2x - p1x;
        dY = p2y - p1y;
        dZ = p2z - p1z;

	temp = dX*dX + dY*dY + dZ*dZ;
        invL = 1.0 / sqrt(temp);

        vec[0] = dX * invL;
        vec[1] = dY * invL;
        vec[2] = dZ * invL;


/*
 *      Build the basis (vec,z,y)
 */	

        vecMinIndex = (fabs(vec[1]) < fabs(vec[0]));
        vecMinIndex = (fabs(vec[2]) < fabs(vec[vecMinIndex])) ? 2 : vecMinIndex;

        xHat[0] = 0.0;
        xHat[1] = 0.0;
        xHat[2] = 0.0;

        xHat[vecMinIndex] = 1.0;

        tmp1 = DotProduct(vec, xHat);
        
        for (i = 0; i < 3; i++) {
            xHat[i] -= tmp1 * vec[i];
        }

        NormalizeVec(xHat);
        cross(vec, xHat, yHat);

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    dGdx[i][j][k] = 0.0;
                }
            }
        }

        for (i = 0; i < nMax; i++) {
            real8 psi, sinPsi, cosPsi;
            real8 z[3];
            real8 mZ[3][3];
            real8 fZ[3][3];
            real8 mStar[3][3];

            psi = M_PI * (i+1) / nMax;

            sinPsi = sin(psi);
            cosPsi = cos(psi);

            for (j = 0; j < 3; j++) {
                z[j] = cosPsi * xHat[j] + sinPsi * yHat[j];
            }

/*
 *          mStar = (zz)^({-1}
 */
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    mZ[j][k] = 0.0;
                }
            }

            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        for (m = 0; m < 3; m++) {
                            mZ[l][k] += C[l][j][k][m] * z[j] * z[m];
                        }
                    }
                }
            }  /* end for (j = 0; j < 3; j++) */

            if (!Matrix33Invert(mZ, mStar)) {
                printf("ERROR: Unabled to invert mZ\n");
                exit(0);
            }

            for (j = 0; j < 3; j++) {
                fZ[j][0] = 0.0;
                fZ[j][1] = 0.0;
                fZ[j][2] = 0.0;
            }

            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        for (m = 0; m < 3; m++) {
                            for (n = 0; n < 3; n++) {
                                for (p = 0; p < 3; p++) {
                                    fZ[j][k] += C[l][m][n][p] *
                                                mStar[j][l] *
                                                mStar[n][k] *
                                                (z[m]*vec[p] + z[p]*vec[m]);
                                }
                            }
                        }
                    }
                }
            }

/*
 *          Integrand
 */
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        dGdx[j][k][l] = dGdx[j][k][l] -
                                        vec[l] * mStar[j][k] +
                                        z[l] * fZ[j][k];
                    }
                }
            }

        }  /* end for (i = 0; i < nMax; i++) */

        tmp1 = M_PI / (real8)nMax;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
		  dGdx[i][j][k] *= tmp1*INV_4_PI_SQUARED*invL*invL;
		}
            }
        }


        for (i = 0; i < 3; i++) 
	  {
	    dGdx2[i][0] =  dGdx[i][0][0];
	    dGdx2[i][1] =  dGdx[i][1][1];
	    dGdx2[i][2] =  dGdx[i][2][2];
	    dGdx2[i][3] = (dGdx[i][1][2] + dGdx[i][2][1]);
	    dGdx2[i][4] = (dGdx[i][2][0] + dGdx[i][0][2]);
	    dGdx2[i][5] = (dGdx[i][1][0] + dGdx[i][0][1]);
	  }

        return;
}



void DerivativeOfGreensFunction2(real8 p1x, real8 p1y, real8 p1z,
				 real8 p2x, real8 p2y, real8 p2z,
				 int qMax,
				 real8 *FqReal[18], real8 *FqImag[18],
				 real8_complex rotMatrix12[3],
				 real8 rotMatrix3[3],
				 real8 dGdx2[3][6])
{
  int   twoqMaxp3 = 2*qMax+3,i,j,k,local,local2;
  int ind, ind1;
  real8 dX, dY, dZ,r4,i4;
  real8 temp, invL, t[3];
  real8 te3[twoqMaxp3], dGdx2aux[18];
  real8_complex te12[twoqMaxp3], J;//[twoqMaxp3];
  real8 coeff;

/*
 *      Define line direction
 */
        dX = p2x - p1x;
        dY = p2y - p1y;
        dZ = p2z - p1z;

	temp = dX*dX + dY*dY + dZ*dZ;
        invL = 1.0 / sqrt(temp);

        t[0] = dX * invL;
        t[1] = dY * invL;
        t[2] = dZ * invL;

	coeff = INV_4_PI_SQUARED*invL*invL;

        te12[0] = 1.0;
        te12[1] = DotProduct(t, rotMatrix12);

        te3[0] = 1.0;
        te3[1] = DotProduct(t, rotMatrix3);


        for (i = 2; i < twoqMaxp3; i++) {
	  te12[i] = te12[i-1] * te12[1];
	  te3[i]  = te3[i-1]  * te3[1];
        }

	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 6; j++) 
	    dGdx2[i][j] = 0.0;

	local=0;
        for (ind = 0; ind < qMax+1; ind++) 
	  for (ind1 = 0; ind1 < 2*ind+2; ind1++) 
	    {
	      J = te12[ind1]*te3[2*ind+1-ind1];
	      for (i = 0; i < 3; i++) 
		for (j = 0; j < 6; j++) 
		  {
		    dGdx2[i][j] += (creal(J)*FqReal[6*i+j][local] - 
				    cimag(J)*FqImag[6*i+j][local])*coeff;
		  }
	      local++;
	    }

}



void DerivativeOfGreensFunction3(real8 p1x, real8 p1y, real8 p1z,
				 real8 p2x, real8 p2y, real8 p2z,
				 int qMax,
				 int *indexList1Real, int *indexList2Real,
				 int *indexList1Imag, int *indexList2Imag,
				 real8 *FqReal, real8 *FqImag,
				 int FqRealNumElem, int FqImagNumElem,
				 real8_complex rotMatrix12[3],
				 real8 rotMatrix3[3],
				 real8 dGdx2[3][6])
{
  int   twoqMaxp3 = 2*qMax+3,i,j,k,local,local2;
  int ind, ind1;
  real8 dX, dY, dZ,r4,i4;
  real8 temp, invL, t[3];
  real8 te3[twoqMaxp3], dGdx2aux[18];
  real8_complex te12[twoqMaxp3], J[(qMax+1)*(qMax+2)];
  real8 coeff;

/*
 *      Define line direction
 */
        dX = p2x - p1x;
        dY = p2y - p1y;
        dZ = p2z - p1z;

	temp = dX*dX + dY*dY + dZ*dZ;
        invL = 1.0 / sqrt(temp);

        t[0] = dX * invL;
        t[1] = dY * invL;
        t[2] = dZ * invL;

	coeff = INV_4_PI_SQUARED*invL*invL;

        te12[0] = 1.0;
        te12[1] = DotProduct(t, rotMatrix12);

        te3[0] = 1.0;
        te3[1] = DotProduct(t, rotMatrix3);


        for (i = 2; i < twoqMaxp3; i++) {
	  te12[i] = te12[i-1] * te12[1];
	  te3[i]  = te3[i-1]  * te3[1];
        }

        for (i = 0; i < 18; i++) dGdx2aux[i] = 0;


	local=0;
        for (ind = 0; ind < qMax+1; ind++) 
	  for (ind1 = 0; ind1 < 2*ind+2; ind1++) 
	    {
	      J[local] = te12[ind1]*te3[2*ind+1-ind1];
	      local++;
	    }



	int offsetJ, offsetG, offsetFq;


        for (offsetFq = 0; offsetFq < FqRealNumElem; offsetFq++) {
            offsetJ = indexList1Real[offsetFq];
            offsetG = indexList2Real[offsetFq];
            dGdx2aux[offsetG] += creal(J[offsetJ])*FqReal[offsetFq];
        }

        for (offsetFq = 0; offsetFq < FqImagNumElem; offsetFq++) {
            offsetJ = indexList1Imag[offsetFq];
            offsetG = indexList2Imag[offsetFq];
            dGdx2aux[offsetG] -= cimag(J[offsetJ])*FqImag[offsetFq];
	}



	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 6; j++) 
	    {
	      local2 = 6*i+j;
	      dGdx2[i][j] = dGdx2aux[local2]*coeff;
	    }
}

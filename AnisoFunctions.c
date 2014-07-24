/* 
 * This file contains functions necessary to run the anisotropic case.
 * It is only a wrapper
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <time.h>

#include "AnisoFunctions.h"


real8 Normal(real8 a[3])
{
        return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
}



void GetCbcc(double c11,double c12, double c44, double C66[6][6], double C[3][3][3][3],
	     double C633[6][3][3])
{
  int i,j,k,l;
  int   i3to6[3][3];


        C66[0][0] = c11;
        C66[0][1] = c12;
        C66[0][2] = c12;
        C66[0][3] = 0.0;
        C66[0][4] = 0.0;
        C66[0][5] = 0.0;

        C66[1][0] = c12;
        C66[1][1] = c11;
        C66[1][2] = c12;
        C66[1][3] = 0.0;
        C66[1][4] = 0.0;
        C66[1][5] = 0.0;

        C66[2][0] = c12;
        C66[2][1] = c12;
        C66[2][2] = c11;
        C66[2][3] = 0.0;
        C66[2][4] = 0.0;
        C66[2][5] = 0.0;

        C66[3][0] = 0.0;
        C66[3][1] = 0.0;
        C66[3][2] = 0.0;
        C66[3][3] = c44;
        C66[3][4] = 0.0;
        C66[3][5] = 0.0;

        C66[4][0] = 0.0;
        C66[4][1] = 0.0;
        C66[4][2] = 0.0;
        C66[4][3] = 0.0;
        C66[4][4] = c44;
        C66[4][5] = 0.0;

        C66[5][0] = 0.0;
        C66[5][1] = 0.0;
        C66[5][2] = 0.0;
        C66[5][3] = 0.0;
        C66[5][4] = 0.0;
        C66[5][5] = c44;

        i3to6[0][0] = 0;
        i3to6[0][1] = 3;
        i3to6[0][2] = 4;

        i3to6[1][0] = 3;
        i3to6[1][1] = 1;
        i3to6[1][2] = 5;

        i3to6[2][0] = 4;
        i3to6[2][1] = 5;
        i3to6[2][2] = 2;


        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        C[i][j][k][l] = C66[i3to6[i][j]][i3to6[k][l]];
                    }
                }
            }
        }

        for (i = 0; i < 3; i++) 
            for (j = 0; j < 3; j++) 
	  {
	    C633[0][i][j] = C[0][0][i][j];
	    C633[1][i][j] = C[1][1][i][j];
	    C633[2][i][j] = C[2][2][i][j];
	    C633[3][i][j] = C[1][2][i][j];
	    C633[4][i][j] = C[2][0][i][j];
	    C633[5][i][j] = C[1][0][i][j];
	  }

}

/*---------------------------------------------------------------------------
 *
 *      Function:    MatrixDblMultComplex
 *      Description: Multiplies a matrix of double values by a matrix of
 *                   complex values.
 *
 *          NOTE: The physical memory layout of the input and output
 *                matrices may be larger than the actual portions
 *                of the matrices being multiplied.  The assumption
 *                is made that the components of matrix <a> involved
 *                in the operation reside in rows zero thru <aRows>-1 and
 *                columns zero thru <aCols>-1, and the corresponding
 *                values for matrices <b> and <c>.
 *
 *          a         Pointer to row 0 column 0 of first matrix to be
 *                    multiplied
 *          aRows     Row order of matrix <a>
 *          aCols     Column order of matrix <a>
 *          aLD       Leading dimension of matrix <a>
 *          b         Pointer to row 0 column 0 of second matrix to be
 *                    multiplied
 *          bCols     Column order of matrix <b>
 *          bLD       Leading dimension of matrix <b>
 *          c         Pointer to location in which to store the
 *                    results of the matrix multiply.  Matrix <c> is
 *                    assumed to be of at last dimensions <aRows> X
 *                    <bCols>.
 *          cLD       Leading dimension of matrix <c>
 *
 *-------------------------------------------------------------------------*/
void MatrixDblMultComplex(double *a, int aRows, int aCols, int aLD,
                          real8_complex *b, int bCols, int bLD,
                          real8_complex *c, int cLD)
{
        int  k, m, n;
        int  aCol, bCol, cCol;
        int  aRow, bRow, cRow;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aRows; m++) {

            aRow = m;
            cRow = m;

            for (n = 0; n < bCols; n++) {

                bCol = n;
                cCol = n;

                cIndex = cRow * cLD + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aCols; k++) {

                    aCol = k;
                    bRow = k;

                    aIndex = aRow * aLD + aCol;
                    bIndex = bRow * bLD + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }

        return;
}


void NormalizeVec(double vec[3])
{
        double a2, a;

        a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        if (a2 > 0.0) {
            a = sqrt(a2);
            vec[0] /= a;
            vec[1] /= a;
            vec[2] /= a;
        }

        return;
}

void cross(double a[3], double b[3], double c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

int Matrix33Invert(double a[3][3], double b[3][3])
{
        int    i, j, k;
        double  p, fmax, fval, eps = 1.0e-20;
        double  tmp[3][3];

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                b[i][j] = (double)(i == j);
                tmp[i][j] = a[i][j];
            }
        }

        for (i = 0; i < 3; i++) {
            fmax = fabs(tmp[i][i]);
            for (j = i+1; j < 3; j++) {
                if (fabs(tmp[j][i]) > fmax) {
                    fmax = fabs(tmp[j][i]);
                    for (k = 0; k < 3; k++) {
                        p = tmp[i][k];
                        tmp[i][k] = tmp[j][k];
                        tmp[j][k] = p;
                        p = b[i][k];
                        b[i][k] = b[j][k];
                        b[j][k] = p;
                    }
                }
            }

            if (fmax < eps) {
#if 0
                printf("Matrix33Invert: fmax < eps, cannot invert!\n");
#endif
                return(0);
            }

            fval = 1.0 / tmp[i][i];

            for (j = 0; j < 3; j++)   {
                tmp[i][j] *= fval;
                b[i][j] *= fval;
            }

            for (k = 0; k < 3; k++) {
                if (k != i) {
                    fval = tmp[k][i];
                    for (j = 0;  j < 3;  j++) {
                        tmp[k][j] -= fval*tmp[i][j];
                        b[k][j] -= fval*b[i][j];
                    }
                }
            }
        }

        return(1);
}

void gqwp(int n, double *xx, double *ww)
{
        int   ti, m, j, k;
        double e1, nn, t, xo, pkm1, pk, t1, pkp1;
        double den, d1, dpn, d2pn, d3pn, d4pn, u, v, h, p, dp, bp, fx, wf;

        m = (n+1)/2;
        e1 = n*(n+1);
        nn = 1.0 - (1.0 - 1.0/n)/(8*n*n);

        for (ti = 0; ti<m; ti++) {
            t = M_PI/(4*n+2)*(4*ti+3);
            xo = nn*cos(t);
            for (j = 0; j<2; j++) {
                pkm1 = 1.0;
                pk = xo;
                for (k = 2; k <= n; k++) {
                    t1 = xo*pk;
                    pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1;
                    pkm1 = pk;
                    pk = pkp1;
                }
                den = 1.0 - xo*xo;
                d1 = n*(pkm1 - xo*pk);
                dpn = d1/den;
                d2pn = (2*xo*dpn - e1*pk)/den;
                d3pn = (4*xo*d2pn + (2.0-e1)*dpn)/den;
                d4pn = (6*xo*d3pn + (6.0-e1)*d2pn)/den;
                u = pk/dpn;
                v = d2pn/dpn;
                h = -u*(1.0 + (u/2)*(v + u*(v*v - u*d3pn/(3*dpn))));
                p = pk + h*(dpn + (h/2)*(d2pn + (h/3)*(d3pn + (h/4)*d4pn)));
                dp = dpn + h*(d2pn + (h/2)*(d3pn + h*d4pn/3));
                h = h - p/dp; xo = xo + h;
            }
            bp = -xo - h;
            fx = d1 - h*e1*(pk+(h/2)*(dpn+(h/3)*(d2pn + (h/4)*(d3pn + (h/5)*d4pn))));
            wf = 2*(1.0 - bp*bp)/(fx*fx);

            xx[ti] = -bp;
            xx[n-1-ti] = bp;
            ww[ti] = wf;
            ww[n-1-ti] = wf;
        }

        if (n%2 > 0) xx[n/2] = 0.0;

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    PreFac
 *      Description: 
 *
 *      Parameters:
 *          IN:     qMax
 *          IN:     binFactTblPtr Pointer to the
 *                          binFacTbl[2*qMax+2][qMax+1][2*qMax+2] array
 *          IN:     glmPtr  Pointer to glm[[2*qMax+2][qMax+2][3][6] array
 *          IN/OUT: FqPtr   Pointer to the Fq[2*qMax+2][qMax+1][3][6] array
 *
 *-------------------------------------------------------------------------*/
void PreFac(int qMax, real8 *binFactTblPtr, real8_complex *glmPtr,
            real8 *ReFqIn, real8 *ImFqIn, 
	    real8 *maxReF, real8 *maxImF)
{
        int i, j, k, l, m, n;
	int ind, ind1,local;
        real8 (*binomialFactCoeffs)[2*qMax+2][2*qMax+2];
        real8_complex (*glm)[2*qMax+2][3][6];
	real8 fac;

        glm = (real8_complex (*)[2*qMax+2][3][6])glmPtr;
        binomialFactCoeffs = (real8 (*)[2*qMax+2][2*qMax+2])binFactTblPtr;

        real8 (*ReFq)[3][(qMax+2)*(qMax+1)];
        ReFq = (real8 (*)[3][(qMax+2)*(qMax+1)])ReFqIn;

        real8 (*ImFq)[3][(qMax+2)*(qMax+1)];
        ImFq = (real8 (*)[3][(qMax+2)*(qMax+1)])ImFqIn;

        for (m = 0; m < 3; m++) 
	  for (n = 0; n < 6; n++) 
	    for (i = 0; i < (qMax+1)*(qMax+2); i++) 
	      {
		ReFq[n][m][i] = 0.0;
		ImFq[n][m][i] = 0.0;
	      }


        for (m = 0; m < 3; m++) {
            for (n = 0; n < 6; n++) {
                for (i = 0; i < qMax+1; i++) {

                    l = 2 * i + 1;

                    //for (j = 0; j < l+1; j++) {
                    for (j = 0; j < 2*i+2; j++) {
                        int intlmj;

			if (j==0)
			  fac = 1;
			else
			  fac = 2;


                        if ((l-j) % 2 == 0) {
                            intlmj = (l-j) / 2;
                        } else {
                            intlmj = (l-j-1) / 2;
                        }

                        for (k = 0; k < intlmj+1; k++) {

			  ind1 = j;
			  ind = i-k;
			  local = ind*(ind+1) + ind1;	    

			  ReFq[n][m][local] += creal(fac * 
						   binomialFactCoeffs[j][l][k] *
						   glm[l][j][m][n]);

			  ImFq[n][m][local] += cimag(fac * 
						   binomialFactCoeffs[j][l][k] *
						   glm[l][j][m][n]);
                        }
                    }
                }
            }
        }


  *maxReF = 0.0;
  *maxImF = 0.0;
  for (m = 0; m < 3; m++) 
    for (n = 0; n < 6; n++) 
      for (i = 0; i < (qMax+2)*(qMax+1); i++) 
	{
	  if (fabs(ReFq[n][m][i]) > *maxReF) *maxReF = fabs(ReFq[n][m][i]);
	  if (fabs(ImFq[n][m][i]) > *maxImF) *maxImF = fabs(ImFq[n][m][i]);
	}

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    GetGreensFuncDerivatives
 *      Description:
 *
 *      Parameters:
 *          IN:     nMax
 *          IN:     C
 *          OUT:    dGdx
 *
 *-------------------------------------------------------------------------*/
void GetGreensFuncDerivatives(int nMax, real8 vec[3], real8 C[3][3][3][3],
                              real8 dGdx[3][3][3])
{
        int   i, j, k, l, m, n, p;
        int   vecMinIndex;
        real8 tmp1;
        real8 xHat[3], yHat[3];
/*
 *      Build the basis (vec,z,y)
 */
        NormalizeVec(vec);

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
                    dGdx[i][j][k] *= tmp1;
                }
            }
        }

        return;
}


void GetGLMandbinom(int qMax, int numTheta, int numPhi, real8 *phi, 
		    real8 *dGdxgridIn,
		    real8 *xx, real8 *ww, 
		    real8_complex glm[2*qMax+2][2*qMax+2][3][6],
		    real8 binomialFactCoeffs[2*qMax+2][2*qMax+2][2*qMax+2])
{
        int   a, b;
        int   i, j, k, l, m, n;
        int   lMax = 2*qMax+1;
        real8 nFac;
        real8 twoPIovernumPhi;
        real8_complex glmaux[lMax+1][lMax+1][3][3][3];
        real8_complex G[numTheta][numPhi];
        real8_complex Ylm[numTheta][numPhi];
        real8_complex Ylmbar[numTheta][numPhi];
        //real8 PlmCol[numTheta];
        real8_complex expPhi[numPhi];
        real8_complex sumYlmbarG[numTheta];
        real8 plm[numTheta][lMax+1];
        real8 plmCol[numTheta];
        real8 q2[lMax+1][lMax+1];
        real8 (*dGdxgrid)[numPhi][3][3][3];


        dGdxgrid = (real8 (*)[numPhi][3][3][3])dGdxgridIn;

        twoPIovernumPhi = 2 * M_PI / numPhi;

	int L,M,intLmM;
	real8 factor;

	binomialFactCoeffs[0][0][0] = 1.0/sqrt(4*M_PI);
	
	
	for (L = 1; L <= lMax; L++) {
	  factor = sqrt((2.0*L+1.0)*(2.0*L-1.0))/(L*1.0);    
	  binomialFactCoeffs[0][L][0] = binomialFactCoeffs[0][L-1][0]*factor;                            
	  
	  for (k = 1; k <= L; k++) {
	    factor = -(L-2.0*k+2.0)*(L-2.0*k+1.0)/(2.0*k*(2.0*L-2.0*k+1.0)); 
	    binomialFactCoeffs[0][L][k] = binomialFactCoeffs[0][L][k-1]*factor;
	  }
	  
	  for (M = 1; M <= L; M++) {
	    factor = -sqrt((2.0*L-1.0)*(2.0*L+1.0)/(1.0*(L+M)*(L+M-1.0)));    
	    binomialFactCoeffs[M][L][0] = binomialFactCoeffs[M-1][L-1][0]*factor;
	    
	    if ((L-M) % 2 == 0)
	      intLmM = (L-M)/2;
	    else
	      intLmM = (L-M-1)/2;
	    
	    for (k = 1; k <= intLmM; k++) {
	      factor = sqrt((2.0*L+1.0)/(2.0*L-1.0))/sqrt(1.0*(L+M)*(L+M-1.0) )* (L-2.0*k-M+2.0)*(L-2.0*k-M+1.0)/(k*2.0);
	      binomialFactCoeffs[M][L][k] = binomialFactCoeffs[M-1][L-1][k-1]*factor;
	    }
	  }
	}
	

         for (l = 0; l <= lMax; l++) 
            for (j = 0; j <= lMax; j++) 
	      q2[l][j] = 0.0;
	   



       for (i = 0; i <= qMax; i++) {

            l = 2*i+1;

            q2[l][0] = 1.0;

            for (j = 1; j <= l; j++) {

                q2[l][j] = q2[l][j-1] / sqrt((real8)(l-j+1) * (real8)(l+j));
            }
        }


        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {

/*
 *                  Copy portion of dGdxgrid into G.  Initially
 *                  only the 'real' portion of G values are set.
 *                  further on the imaginary portions may be updated.
 */
                    for (m = 0; m < numTheta; m++) {
                        for (n = 0; n < numPhi; n++) {
                            G[m][n] = dGdxgrid[m][n][i][j][k];
                        }

                    }

                    for (m = 0; m <= lMax; m++) {
                        real8 nFac0;
/*
 *                      Plm is a numTheta X m+1 matrix (to max of
 *                      numTheta X lMax+1)
 */
                        for (n = 0; n < numTheta; n++) {
                            ComputeLegendre(m, xx[n], plm[n]);
                        }

                        nFac0 = sqrt((2.0 * (real8)m + 1.0) / fourPI);
                        for (n = 0; n <= m; n++) {

			  nFac = nFac0 * q2[m][abs(n)];

                            for (a = 0; a < numTheta; a++) {
                                plmCol[a] = plm[a][abs(n)];
                            }
                            for (b = 0; b < numPhi; b++) {
                                expPhi[b] = cos(n*phi[b]) + sin(n*phi[b]) * I;
                            }

                            MatrixDblMultComplex(plmCol, numTheta, 1, 1,
                                                 expPhi, numPhi, numPhi,
                                                 (real8_complex*)Ylm, numPhi);

                            for (a = 0; a < numTheta; a++) {
                                for (b = 0; b < numPhi; b++) {
                                    Ylm[a][b] *= nFac;
                                    Ylmbar[a][b] = creal(Ylm[a][b]) -
                                                   cimag(Ylm[a][b])*I;
                                }
                            }

                            glmaux[m][n][i][j][k] = 0.0;

                            for (a = 0; a < numTheta; a++) {
                                sumYlmbarG[a] = 0.0;
                                for (b = 0; b < numPhi; b++) {
                                    sumYlmbarG[a] += Ylmbar[a][b] * G[a][b] *
                                                     twoPIovernumPhi;
                                }

                                glmaux[m][n][i][j][k] += sumYlmbarG[a] * ww[a];

                            }
#if 0
			    printf("glmaux[%d][%d][%d][%d][%d] = %.15e %.15e\n", m+1, n+1, i+1, j+1, k+1,
				    creal(glmaux[m][n][i][j][k]), cimag(glmaux[m][n][i][j][k]));
#endif
			        for (a = 0; a < numTheta; a++) {
			      for (b = 0; b < numPhi; b++) {
				G[a][b] -= glmaux[m][n][i][j][k] *
				  Ylm[a][b];
			      }
			      }

                        }  /* end for n ... */
                    }  /* end for m ... */
                }  /* end for k ... */
            }  /* end for j ... */
        }  /* end for i ... */

        for (i = 0; i < 3; i++) {
            for (j = 0; j <= lMax; j++) {
                for (k = 0; k <= j; k++) {
                    glm[j][k][i][0] = glmaux[j][k][i][0][0];
                    glm[j][k][i][1] = glmaux[j][k][i][1][1];
                    glm[j][k][i][2] = glmaux[j][k][i][2][2];
                    glm[j][k][i][3] = glmaux[j][k][i][1][2] +
                                      glmaux[j][k][i][2][1];
                    glm[j][k][i][4] = glmaux[j][k][i][2][0] +
                                      glmaux[j][k][i][0][2];
                    glm[j][k][i][5] = glmaux[j][k][i][0][1] +
                                      glmaux[j][k][i][1][0];
#if 0
		    printf("glm[%d][%d][%d][*]=%6.3f updated\n", j+1, k+1, i+1,glm[j][k][i][0]);
#endif
                }
            }
        }
	
        return;
}


void CreateLists(int qMax, real8 maxF,real8 *FqIn,
		 int *ind1list, int *ind2list,
		 real8 *FqArray, int *sizeFq)
{
  int l,m,i,j;
  int ind;
  int local;

  real8 (*Fq)[3][(qMax+2)*(qMax+1)];
  Fq = (real8 (*)[3][(qMax+2)*(qMax+1)])FqIn;


  ind = 0;
  for (m = 0; m < qMax+1; m++) 
    for (l = 0; l < 2*m+2; l++) 
      {
	local = m*(m+1)+l;
	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 6; j++) 
	    if ( fabs(Fq[j][i][local]) >  maxF*1e-15) 
	      {
		ind1list[ind]=local;
		ind2list[ind]=6*i+j;
		FqArray[ind]=Fq[j][i][local];
		ind++;
	      }
      }

  *sizeFq = ind;
}

#if 0
void CreateLists3(int qMax, real8 maxF,real8 *FqIn,real8 *FqArray[18])
{
  int l,m,i,j;
  int ind;
  int local;

  real8 (*Fq)[3][(qMax+2)*(qMax+1)];
  Fq = (real8 (*)[3][(qMax+2)*(qMax+1)])FqIn;

  for (i = 0; i < 3; i++) 
    for (j = 0; j < 6; j++) 
      {
	for (local = 0; local < (qMax+2)*(qMax+1); local++) 
	  {
	    FqArray[6*i+j][local]=Fq[j][i][local];
	  }
      }
}
#endif


/* To access the P_l,m in a array */
/* {0,0}{1,0}{1,1}{2,0}{2,1}{2,2} ... */
int atLm(int l, int m){
    /* summation series over l + m => (l*(l+1))/2 + m */
    return ((l*(l+1))>>1) + m;
}
/*
 *      For use when calculating the Legendre polynomials; to
 *      access P_l,m in an array.
 *
 *      {0,0}{1,0}{1,1}{2,0}{2,1}{2,2} ...
 */
#define ATLM(l,m) ( (((l)*((l)+1))>>1) + (m) )


/*---------------------------------------------------------------------------
 *
 *      Function:    ComputeLegendre
 *      Description: Compute the legendre polynomial by recurrence.
 *
 *      Parameters:
 *          IN:  P         Degree of the polynomial
 *          IN:  x         Value for which the polynomial is being calculated
 *          OUT: legendre  Array of length P+1 in which to return the
 *                         legendre function for <x>.
 *
 *-------------------------------------------------------------------------*/
void ComputeLegendre(int P, real8 x, real8 *legendre){
        int   l, m;
        real8 factor = -sqrt(1.0-pow(x,2));
        real8 lTemp[((P+2)*(P+1))/2];

/*
 *      Init legendre
 */
        lTemp[atLm(0,0)] = 1.0;        /* P_0,0(x) = 1 */
        lTemp[atLm(1,0)] = x;          /* P_1,0(x) = x */
        lTemp[atLm(1,1)] = factor;     /* P_1,1(x) = -sqrt(1 - x^2) */

        for (l = 2; l <= P ; ++l ){
            for (m = 0; m < l - 1 ; ++m ){
                /* P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k) */
	      lTemp[atLm(l,m)] = ((real8)(2*l-1) * x * lTemp[atLm(l-1,m)] -
                                    (real8)(l + m - 1) * lTemp[atLm(l-2,m)]) /
                                   (real8)(l-m);
            }

            /* P_l,l-1 = (2l-1)*x*P_l-1,l-1 */
            lTemp[atLm(l,l-1)] = (real8)(2*l-1) * x * lTemp[atLm(l-1,l-1)];

            /* P_l,l = (2l-1)*factor*P_l-1,l-1 */
            lTemp[atLm(l,l)] = (real8)(2*l-1) * factor * lTemp[atLm(l-1,l-1)];
        }

        for (m = 0 ; m <= P ; ++m){
            legendre[m] = lTemp[atLm(P,m)];
        }

        return;
}



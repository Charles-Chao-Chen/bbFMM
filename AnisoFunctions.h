#ifndef _ANISOFUNCTIONS
#define _ANISOFUNCTIONS

#define DotProduct(vec1,vec2)	    \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))

#include <complex.h>

typedef double real8;
typedef complex double real8_complex;

#define fourPI 1.256637061435917e+01
#define INV_4_PI_SQUARED 0.025330295910584  /* 1/(4*pi^2) */


real8 Normal(real8 a[3]);
void GetCbcc(double c11, double c12, double c44, double C66[6][6], 
	     double C[3][3][3][3], double C633[6][3][3]);
void MatrixDblMultComplex(double *a, int aRows, int aCols, int aLD,
                          real8_complex *b, int bCols, int bLD,
                          real8_complex *c, int cLD);
void NormalizeVec(double vec[3]);
  void cross(double a[3], double b[3], double c[3]);
  int Matrix33Invert(double a[3][3], double b[3][3]);

void gqwp(int n, double *xx, double *ww);

void PreFac(int qMax, double *binFactTblPtr, real8_complex *glmPtr,
            double *ReFqIn, double *ImFqIn, 
	    double *maxReF, double *maxImF);


void GetGreensFuncDerivatives(int nMax, real8 vec[3], real8 C[3][3][3][3],
                              real8 dGdx[3][3][3]);


void DerivativeOfGreensFunction(int nMax, real8 C[3][3][3][3],
				real8 p1x, real8 p1y, real8 p1z,
				real8 p2x, real8 p2y, real8 p2z,
				real8 dGdx[3][6]);

void DerivativeOfGreensFunction2(real8 p1x, real8 p1y, real8 p1z,
				 real8 p2x, real8 p2y, real8 p2z,
				 int qMax,
				 real8 *FqReal[18], real8 *FqImag[18],
				 real8_complex rotMatrix12[3],
				 real8 rotMatrix3[3],
				 real8 dGdx2[3][6]);

void DerivativeOfGreensFunction3(real8 p1x, real8 p1y, real8 p1z,
				 real8 p2x, real8 p2y, real8 p2z,
				 int qMax,
				 int *indexList1Real, int *indexList2Real,
				 int *indexList1Imag, int *indexList2Imag,
				 real8 *FqReal, real8 *FqImag,
				 int FqRealNumElem, int FqImagNumElem,
				 real8_complex rotMatrix12[3],
				 real8 rotMatrix3[3],
				 real8 dGdx2[3][6]);


void GetGLMandbinom(int qMax, int numTheta, int numPhi, real8 *phi, 
		    real8 *dGdxgridIn,
		    real8 *xx, real8 *ww, 
		    real8_complex glm[2*qMax+2][2*qMax+2][3][6],
		    real8 binomialFactCoeffs[2*qMax+2][2*qMax+2][2*qMax+2]);

void CreateLists(int qMax, real8 maxF,real8 *FqIn,
		 int *ind1list, int *ind2list,
		 real8 *FqArray, int *sizeFq);

void ComputeLegendre(int P, real8 x, real8 *legendre);



void IntegrandStressAtPoint(real8 x, real8 y, real8 z,
			    real8 px, real8 py, real8 pz,      
			    real8 t[3],
			    real8 bx,  real8 by,  real8 bz,
			    real8 a, int qMax,
			    real8 C[6][6],
			    int *indexList1Real, int *indexList2Real,
			    int *indexList1Imag, int *indexList2Imag,
			    real8 *FqReal, real8 *FqImag,
			    int FqRealNumElem, int FqImagNumElem,
			    real8_complex rotMatrix12[3],
			    real8 rotMatrix3[3],
			    real8 dGdx[18],real8 dsigma[6]);




void CheckIntegrand(real8 pointx, real8 pointy, real8 pointz,
		    real8 p1x, real8 p1y, real8 p1z,
		    real8 p2x, real8 p2y, real8 p2z,
		    real8 bx,  real8 by,  real8 bz,
		    real8 a,
		    int qMax,real8 C66[6][6],
		    int *ind1Re, int *ind2Re,
		    int *ind1Im, int *ind2Im,
		    real8 *ReFqout, real8 *ImFqout,
		    int sizeReFq, int sizeImFq,
		    real8_complex rotMatrix12[3],
		    real8 rotMatrixRow3[3],
		    real8 sigma[6]);

void ChaoStress(double p1x, double p1y, double p1z,
		double p2x, double p2y, double p2z,      
		int qMax,
		int *indexList1Real, int *indexList2Real,
		int *indexList1Imag, int *indexList2Imag,
		real8 *FqReal, real8 *FqImag,
		int FqRealNumElem, int FqImagNumElem,
		real8_complex e12[3],real8 e3[3], real8 C633[6][3][3],
                double TrapdGdx[3][6],double K[54]);
  
#endif

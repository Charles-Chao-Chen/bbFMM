#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "directCalc.h"


void EvaluateLayer(vec3 *field, int Nf, vec3 *source, double *intensity,
		   int Ns, int2 dof, int lpbc, double boxLen,
		   kfun_t kfun, double *potential) {
   
  int i, j, l1, l2, l3;

  // for dgemv
  int incr=1;
  char trans[] = "n";
  double alpha = 1, beta = 1;
	      
  int start = floor(0.5 + (pow(3,lpbc)-1)/2.);      
  int end;
  if (lpbc <  2)
    end  = -1;
  else
    end  =  1;

  // initialize the potential
  int dofNf = dof.f*Nf;
  for (i=0; i<dofNf; i++)
    potential[i] = 0;
  
  // Compute the interactions from 'start' level to 'end' level
  vec3 sourcepos, cshift;
  double Kij[dof.f*dof.s];

  for (i=0;i<Nf;i++) {
    for (j=0;j<Ns;j++) {
      
      for (l1=-start;l1<start+1;l1++) {
	cshift.x = (double)l1*boxLen;
	for (l2=-start;l2<start+1;l2++) {
	  cshift.y = (double)l2*boxLen;
	  for (l3=-start;l3<start+1;l3++) {
	    cshift.z = (double)l3*boxLen;

	    // apply the shift of source points
	    if (abs(l1) > end || abs(l2) > end || abs(l3) > end) {
	      sourcepos.x = source[j].x + cshift.x;
	      sourcepos.y = source[j].y + cshift.y;
	      sourcepos.z = source[j].z + cshift.z;
              
	      kfun(&field[i], &sourcepos, Kij);

	      if (dof.f == 1 && dof.s == 1)
		potential[i] += Kij[0] * intensity[j];
	      else
		dgemv_(trans, &dof.f, &dof.s, &alpha, Kij, &dof.f, intensity + j*dof.s, &incr,
		       &beta, potential + i*dof.f, &incr);
	    }
	  }
	}
      }
    }
  }
}


/*
  // Old function: explicitly creates a big coefficient matrix and use dgemv
  // which can result in memory usage problem
void EvaluateLayer(vec3 *field, int Nf, vec3 *source, double *sourceCharge,
		   int Ns, int2 dof, int lpbc, double boxLen,
		   kfun_t kfun, double *stressField) {
   
  int i, j, k, l, count, l1, l2, l3, begGlobal, incr=1;
  int dofNf = dof.f*Nf, dofNs = dof.s*Ns;
   
  int start = floor(0.5 + (pow(3,lpbc)-1)/2.);      
  int end;
  if (lpbc <  2)
    end  = -1;
  else
    end  =  1;

  // Compute the interactions from 'start' level to 'end' level
  vec3 sourcepos, cshift;
  double Kij[dof.f*dof.s];
  //double *KPBC = (double *) calloc(dofNf*dofNs, sizeof(double));
  double *KPBC = (double *) calloc(dofNf*dofNs, sizeof(double));

  for (i=0;i<Nf;i++) {
    for (j=0;j<Ns;j++) {
      
      begGlobal = j*dofNf*dof.s + i*dof.f;
      
      for (l1=-start;l1<start+1;l1++) {
	cshift.x = (double)l1*boxLen;
	for (l2=-start;l2<start+1;l2++) {
	  cshift.y = (double)l2*boxLen;
	  for (l3=-start;l3<start+1;l3++) {
	    cshift.z = (double)l3*boxLen;
	    if (abs(l1) > end || abs(l2) > end || abs(l3) > end) {
	      sourcepos.x = source[j].x + cshift.x;
	      sourcepos.y = source[j].y + cshift.y;
	      sourcepos.z = source[j].z + cshift.z;
              
	      //EvaluateKernel(field[i], sourcepos, kfun, Kij, dof);
	      kfun(&field[i], &sourcepos, Kij);
		  
	      count = 0;
	      for (l=0; l<dof.s; l++)
		for (k=0; k<dof.f; k++, count++)
		  KPBC[begGlobal + l*dofNf + k] += Kij[count];
	    }
	  }
	}
      }
    }
  }

  char trans[] = "n";
  double alpha = 1, beta = 0;
  dgemv_(trans, &dofNf, &dofNs, &alpha, KPBC, &dofNf, sourceCharge, &incr,
	 &beta, stressField, &incr);
  
  free(KPBC), KPBC=NULL; 
}
*/


// Computes mean stress over the box directly
void MeanStressDirect(int numGauss, vec3 *source, double
		      *sourceCharge, int sourceSize, int2 dof, int lpbc,
		      double Len, kfun_t kfun, double
		      *meanStress){
   
   
  double gpoint[numGauss], gweight[numGauss];
  gqwp(numGauss, gpoint, gweight);

  // Initialize mean stress
  int idof;
  for (idof=0; idof<dof.f; idof++)
    meanStress[idof] = 0;
   
  int i, j, k;
  double halfLen = Len/2;
  for (i=0; i<numGauss; i++)
    for (j=0; j<numGauss; j++)
      for (k=0; k<numGauss; k++) {
	vec3 gpoint3D;
	gpoint3D.x = gpoint[i]*halfLen;
	gpoint3D.y = gpoint[j]*halfLen;
	gpoint3D.z = gpoint[k]*halfLen;
         
	double stress[dof.f];
	EvaluateLayer(&gpoint3D, 1, source, sourceCharge, sourceSize,
		      dof, lpbc, Len, kfun, stress);

	// Update mean stress
	double gweight3D = gweight[i]*gweight[j]*gweight[k];
	for (idof=0; idof<dof.f; idof++)
	  meanStress[idof] += gweight3D * stress[idof];
      }

  // divide by the volume of the cubic
  for (idof=0; idof<dof.f; idof++)
    meanStress[idof] /= 8;
}


/*
 * Function: DirectCalc3D
 * ---------------------------------------------------------------------
 * Computes the potential at the first field point and returns 
 * the result in phi.
 */
/*
void DirectCalc3D(vec3 *field, int Nf, vec3 *source, double
		  *sourceCharge, int sourceSize, int2 dof, int levelpbc,
		  double boxLen, kfun_t kfun, double
		  *stressField) {
  */
double* directCalc( vec3 *field, int Nf, vec3 *source, int sourceSize,
		    double *sourceCharge,
		    int2 dof, double boxLen, kfun_t kfun, int levelpbc ) {

  double *stressField = calloc(Nf*dof.f, sizeof(double));
  
  // Read direct result from file
  FILE *pDirFile;
  char filename[50];
  sprintf(filename, "DirectResultlpbc%d.bin", levelpbc);
  pDirFile = fopen(filename, "rb");
  
  if (pDirFile != NULL) {

    int i = fread(stressField, sizeof(double), Nf*dof.f, pDirFile);
    fclose(pDirFile);

    if (i != Nf*dof.f)
      printf("fread() error in DirectCalc3D().\n");
    else
      printf("reading direct result from file.\n");
  }
  else {
    //fclose(pDirFile);	
	
    if (levelpbc < 2)
      EvaluateLayer(field, Nf, source, sourceCharge, sourceSize, dof,
		    levelpbc, boxLen, kfun, stressField);
    
        
    else { // levelpbc >= 2
      
      // compute stress from 3x3 domain
      EvaluateLayer(field, Nf, source, sourceCharge, sourceSize, dof, 1,
		    boxLen, kfun, stressField);
      
      // compute PBC
      int i;
      double *stressFieldPBC;
      stressFieldPBC = malloc( Nf*dof.f * sizeof(double));
      
      EvaluateLayer(field, Nf, source, sourceCharge, sourceSize, dof,
		    levelpbc, boxLen, kfun, stressFieldPBC);
      
      
      // Calculate mean stress using Gauss quadrature
      int numGauss = 10;
      double meanStressDir[dof.f];
      MeanStressDirect(numGauss, source, sourceCharge, sourceSize,
		       dof, levelpbc, boxLen, kfun,
		       meanStressDir);
      
      
      // Subtract mean stress from pbc stress
      int j;
      for (i=0; i<Nf; i++)
	for (j=0; j<dof.f; j++)
	  stressFieldPBC[i*dof.f + j] -= meanStressDir[j];
      

      print_array(stressField, Nf * dof.f, "3x3x3 domain");
      print_array(stressFieldPBC, Nf * dof.f, "Direct PBC");
	  
      // Update field stess
      int fieldStressSize = Nf * dof.f;
      for (i=0; i < fieldStressSize; i++)
	stressField[i] += stressFieldPBC[i];
      //stressField[i] = stressFieldPBC[i];
      
      free(stressFieldPBC), stressFieldPBC = NULL;
    }


    /*
    // Write result to file
    //FILE *pDirFile;
    //char filename[50];
    //sprintf(filename, "DirectResultlpbc%d.bin", levelpbc);
    pDirFile = fopen(filename, "wb");
  
    if (pDirFile != NULL) {
      fwrite(stressField, sizeof(double), Nf*dof.f, pDirFile);
      fclose(pDirFile);
      printf("writing direct result to file.\n");
    }
    else
      printf("Can't open file to write direct result!\n");
    */
  }

  return stressField;
}


/*
 * Function: EwaldSolution
 * -------------------------------------------------------------------
 * Computes the Ewald solution for 1/r, r.x/r^3, and 1/r^4 kernels
 * (for comparision to PBC with FMM).
 */
void EwaldSolution(vec3 *field, vec3 *source, double *q, int Nf, 
		   int Ns, int nstart, int mstart, double L, double beta,
		   double *phi, double *corr) {
  /*******************************
   *    LAPLACIAN KERNEL (1/r)   *
   *******************************/
#ifdef LAPLACIAN
  int i, j, m1, m2, m3, n1, n2, n3;
  vec3 diff, cshift;
  double sum, tot, r, s, s2, f, g, mdotr;
  double pi = M_PI;
  double beta2 = beta*beta;
	
  tot = 0;
  for (i=0;i<Nf;i++) {
    sum = 0;
    for (j=0;j<Ns;j++) {
      // Compute the direct contribution for n = [0 0 0]
      if (source[j].x != field[i].x || source[j].y != field[i].y ||
	  source[j].z != field[i].z) {
	diff.x = source[j].x - field[i].x;
	diff.y = source[j].y - field[i].y;
	diff.z = source[j].z - field[i].z;
	r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
				
	s = beta*r;
	sum += q[j]*erfc(s)/r;
      }
			
      // Compute the rest of the direct sum
      for (n1=-nstart;n1<nstart+1;n1++) {
	cshift.x = (double)n1*L;
	for (n2=-nstart;n2<nstart+1;n2++) {
	  cshift.y = (double)n2*L;
	  for (n3=-nstart;n3<nstart+1;n3++) {
	    cshift.z = (double)n3*L;
	    if (n1 != 0 || n2 != 0 || n3 != 0) {
	      diff.x = (source[j].x + cshift.x) - field[i].x;
	      diff.y = (source[j].y + cshift.y) - field[i].y;
	      diff.z = (source[j].z + cshift.z) - field[i].z;
	      r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
							
	      s = beta*r;
	      sum += q[j]*erfc(s)/r;
	    }
	  }
	}
      }
    }
		
    // Add the direct contribution and subtract self-energy to the total energy
    tot += q[i]*(sum - 2.0*beta*q[i]/sqrt(pi));
  }
	
  sum = 0;
  // Compute the reciprocal sum
  for (m1=-mstart;m1<mstart+1;m1++) {
    for (m2=-mstart;m2<mstart+1;m2++) {
      for (m3=-mstart;m3<mstart+1;m3++) {
	if (m1 != 0 || m2 != 0 || m3 != 0) {
	  r = sqrt((double)(m1*m1+m2*m2+m3*m3));	      
	  s = pi*r/beta;
	  s2 = s*s;
	  f = 0;
	  g = 0;
	  for (j=0;j<Ns;j++) {
	    mdotr = (double)m1*source[j].x + (double)m2*source[j].y +
	      (double)m3*source[j].z;
	    f += q[j]*cos(2*pi*mdotr);
	    g += q[j]*sin(2*pi*mdotr);
	  }
	  sum += 1.0/s2*(f*f+g*g)*exp(-s2);
	}
      }
    }
  }
	
  // Add the reciprocal contribution to the total energy
  tot += sum*pi/beta2;
	
  // Return potential energy
  *phi = tot;
	
  // Compute the spherical shell extrinisic correction
  diff.x = 0;
  diff.y = 0;
  diff.z = 0;
  for (j=0;j<Ns;j++) {
    diff.x += q[j]*source[j].x;
    diff.y += q[j]*source[j].y;
    diff.z += q[j]*source[j].z;
  }
	
  // Return correction to the total energy (for comparision with direct sum)
  *corr = 4.0*pi/3.0*(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
	
	
  /****************************************************
   * LAPLACIAN FORCE (r.x/r^3) - only the x-component *
   ****************************************************/
#elif LAPLACIANFORCE
  int i, j, m1, m2, m3, n1, n2, n3;
  vec3 diff, cshift;
  double sum, r, r2, s, s2, f, g, mdotr;
  double pi = M_PI;
  double fac1 = 2.0*beta/sqrt(pi), fac2;
  double *fi;
  fi = (double *)malloc(Nf * sizeof(double));
	
  for (i=0;i<Nf;i++) {
    sum = 0;
    for (j=0;j<Ns;j++) {
      // Compute force contribution from the n = [0 0 0] term of direct sum
      if (source[j].x != field[i].x || source[j].y != field[i].y ||
	  source[j].z != field[i].z) {
	diff.x = source[j].x - field[i].x;
	diff.y = source[j].y - field[i].y;
	diff.z = source[j].z - field[i].z;
	r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
	r2 = r*r;
	s = beta*r;
	s2 = s*s;
	sum += q[j]*(fac1*exp(-s2)+erfc(s)/r)*diff.x/r2;
      }
			
      // Compute contribution from the rest of the direct sum
      for (n1=-nstart;n1<nstart+1;n1++) {
	cshift.x = (double)n1*L;
	for (n2=-nstart;n2<nstart+1;n2++) {
	  cshift.y = (double)n2*L;
	  for (n3=-nstart;n3<nstart+1;n3++) {
	    cshift.z = (double)n3*L;
	    if (n1 != 0 || n2 != 0 || n3 != 0) {
	      diff.x = (source[j].x + cshift.x) - field[i].x;
	      diff.y = (source[j].y + cshift.y) - field[i].y;
	      diff.z = (source[j].z + cshift.z) - field[i].z;
	      r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
	      r2 = r*r;
	      s = beta*r;
	      s2 = s*s;
	      sum += q[j]*(fac1*exp(-s2)+erfc(s)/r)*diff.x/r2;
	    }
	  }
	}
      }
    }
		
    // Stores the force contribution from the direct sum
    fi[i] = sum;
  }
	
  // Compute the force contribution from the reciprocal sum
  sum = 0;
  for (m1=-mstart;m1<mstart+1;m1++) {
    for (m2=-mstart;m2<mstart+1;m2++) {
      for (m3=-mstart;m3<mstart+1;m3++) {
	if (m1 != 0 || m2 != 0 || m3 != 0) {
	  r = sqrt((double)(m1*m1+m2*m2+m3*m3));
	  r2 = r*r;
	  s = pi*r/beta;
	  s2 = s*s;
	  f = 0;
	  g = 0;
	  for (j=0;j<Ns;j++) {
	    mdotr = (double)m1*source[j].x + (double)m2*source[j].y +
	      (double)m3*source[j].z;
	    f += q[j]*cos(2.0*pi*mdotr);
	    g += q[j]*sin(2.0*pi*mdotr);
	  }
					
	  fac2 = 2.0*(double)m1/r2*exp(-s2);
	  for (i=0;i<Nf;i++) {
	    mdotr = (double)m1*field[i].x + (double)m2*field[i].y +
	      (double)m3*field[i].z;
	    fi[i] += fac2*(g*cos(2.0*pi*mdotr)-f*sin(2.0*pi*mdotr));
	  }
	}
      }
    }
  }
	
  // Return the force vector
  for (i=0;i<Nf;i++)
    phi[i] = q[i]*fi[i];
	
  free(fi);
	
  // Compute the force contribution from spherical shell extrinisic correction
  // First compute the x-component of the dipole
  diff.x = 0;
  for (j=0;j<Ns;j++) {
    diff.x += q[j]*source[j].x;
  }
	
  // Return correction to the force (for comparision with direct sum)
  *corr = 4.0*pi/3.0*diff.x;
	
  /**********************
   *    1/r^4 KERNEL    *
   **********************/
#elif ONEOVERR4
  int i, j, m1, m2, m3, n1, n2, n3;
  vec3 diff, cshift;
  double sum, tot, r, s, s2, s4, f, g, mdotr;
  double pi = M_PI;
  double fac = 2*pi*sqrt(pi)*beta;
  double beta2 = beta*beta;
  double beta4 = beta2*beta2;
	
  tot = 0;
  for (i=0;i<Nf;i++) {
    sum = 0;
    for (j=0;j<Ns;j++) {
      // Compute the direct contribution for n = [0 0 0]
      if (source[j].x != field[i].x || source[j].y != field[i].y ||
	  source[j].z != field[i].z) {
	diff.x = source[j].x - field[i].x;
	diff.y = source[j].y - field[i].y;
	diff.z = source[j].z - field[i].z;
	r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
				
	s = beta*r;
	s2 = s*s;
	s4 = s2*s2;
	g = exp(-s2)*(s2+1.0);
	sum += q[j]*g/s4;
      }
			
      // Compute the rest of the direct sum
      for (n1=-nstart;n1<nstart+1;n1++) {
	cshift.x = (double)n1*L;
	for (n2=-nstart;n2<nstart+1;n2++) {
	  cshift.y = (double)n2*L;
	  for (n3=-nstart;n3<nstart+1;n3++) {
	    cshift.z = (double)n3*L;
	    if (n1 != 0 || n2 != 0 || n3 != 0) {
	      diff.x = (source[j].x + cshift.x) - field[i].x;
	      diff.y = (source[j].y + cshift.y) - field[i].y;
	      diff.z = (source[j].z + cshift.z) - field[i].z;
	      r = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
							
	      s = beta*r;
	      s2 = s*s;
	      s4 = s2*s2;
	      g = exp(-s2)*(s2+1.0);
	      sum += q[j]*g/s4;
	    }
	  }
	}
      }
			
      // Compute the reciprocal sum
      diff.x = source[j].x - field[i].x;
      diff.y = source[j].y - field[i].y;
      diff.z = source[j].z - field[i].z;
      for (m1=-mstart;m1<mstart+1;m1++) {
	for (m2=-mstart;m2<mstart+1;m2++) {
	  for (m3=-mstart;m3<mstart+1;m3++) {
	    r = sqrt((double)(m1*m1+m2*m2+m3*m3));
						
	    s = pi*r/beta;
	    s2 = s*s;
	    f = fac*(exp(-s2)-sqrt(pi)*s*erfc(s));
	    mdotr = (double)m1*diff.x + (double)m2*diff.y + (double)m3*diff.z;
	    sum += q[j]*f*cos(-2*pi*mdotr);
	  }
	}
      }
    }
		
    // Add to total potential energy and subtract self-energy
    tot += q[i]*(sum - 0.5*beta4*q[i]);
  }
	
  // Return the potential energy and correction
  *phi = tot;
  *corr = 0;
	
#endif
} 




void DirWo1 (double *phi, int N, int lpbc) {

  double *a0 = malloc( N*sizeof(double) );
  double *a1 = malloc( N*sizeof(double) );
  double *a2 = malloc( N*sizeof(double) );
  double *a3 = malloc( N*sizeof(double) );
  FILE *f0 = fopen("DirectResultlpbc0.bin", "rb");
  FILE *f1 = fopen("DirectResultlpbc1.bin", "rb");
  FILE *f2 = fopen("DirectResultlpbc2.bin", "rb");
  FILE *f3 = fopen("DirectResultlpbc3.bin", "rb");

  READ_CHECK( fread( a0, sizeof(double), N, f0 ), N );
  READ_CHECK( fread( a1, sizeof(double), N, f1 ), N );
  READ_CHECK( fread( a2, sizeof(double), N, f2 ), N );
  READ_CHECK( fread( a3, sizeof(double), N, f3 ), N );

  fclose(f0);
  fclose(f1);
  fclose(f2);
  fclose(f3);

  phi = Add(phi, a0, N);
  //print_array(a0, N, "PBC 0");
  
  if (lpbc == 2) {
    a2  = Subtract( a2, N, a1, N);
    phi = Add(phi, a2, N);
    //print_array(a2, N, "PBC 2");
    /*
    int i;
    printf("a2 / a0:\n");
    for (i=0; i<N; i++)
      printf("%f\t", a2[i]/a0[i]);
    printf("\n");
    */
  }
  else if (lpbc == 3) {
    a3  = Subtract( a3, N, a1, N);
    phi = Add(phi, a3, N);
  }
  

  free(a0);
  free(a1);
  free(a2);
  free(a3);
}


void DirPBC (double *phi, int N, int lpbc) {

  double *a1 = malloc( N*sizeof(double) );
  double *a2 = malloc( N*sizeof(double) );
  double *a3 = malloc( N*sizeof(double) );
  FILE *f1 = fopen("DirectResultlpbc1.bin", "rb");
  FILE *f2 = fopen("DirectResultlpbc2.bin", "rb");
  FILE *f3 = fopen("DirectResultlpbc3.bin", "rb");

  READ_CHECK( fread( a1, sizeof(double), N, f1 ), N );
  READ_CHECK( fread( a2, sizeof(double), N, f2 ), N );
  READ_CHECK( fread( a3, sizeof(double), N, f3 ), N );

  fclose(f1);
  fclose(f2);
  fclose(f3);

  if (lpbc == 2) {
    a2 = Subtract( a2, N, a1, N);
    phi = Add(phi, a2, N);
  }
  else if (lpbc == 3) {
    a3 = Subtract( a3, N, a1, N);
    phi = Add(phi, a3, N);
  }

  free(a1);
  free(a2);
  free(a3);
}

#include "bbfmm.h"
#include "directCalc.h"
#include <time.h>        // for 'time()'

void run_test(int n, kernel_t kernel, gridT grid) {

  int N          = 1000; // number of particles
  int level_tree = 3;    // hierarchical tree level
  int level_pbc  = 0;    /* lpbc  = 0, 1, 2, 3
			  * start = 0, 1, 4, 13
			  * Note: the original computational domain is lpbc=0
			  *       and PBC starts from lpbc=2
			  */
  double box_len = 4.;
  double alpha   = 0.;

  
  int2 dof = {1, 1};

  // randomly set source and field points
  int Ns = N;
  int Nf = N;
  vec3 *source = malloc( Ns       * sizeof(vec3)   );
  vec3 *field  = malloc( Nf       * sizeof(vec3)   );
  double *q    = malloc( Ns*dof.s * sizeof(double) );

  srand(time(NULL));
  SetSources( field, Nf, source, Ns, q, dof.s, box_len );

  double epsilon = 1e-14;
  /*
  FMMSrc fmm_src = {.srcT   = PNT,
		    .Ns     = Ns,
		    .q      = q,
		    .source = source};
  */
  FMMSrc fmm_src = {PNT, source, NULL, q, NULL, Ns, 0};
  double *phiFmm = bbfmm( fmm_src, field, Nf, dof, box_len, alpha,
			  n, grid, kernel, level_tree, level_pbc,
			  epsilon );

  double *phiDir = directCalc( field, Nf, source, Ns, q, dof, box_len,
			       kernel.kfun, level_pbc );
	
  double err = ComputeL2Err(phiFmm, phiDir, Nf*dof.f);
  

  free(field);
  free(source);
  free(q);
  free(phiFmm);
  free(phiDir);

  printf("Error: %.5f.\n", err);
  assert(err < 1e-14);
}

int main(int argc, char *argv[]) {

  
  DEBUG("Kernel: f=1 (Unif)");
  {
    int n = 1;
    kernel_t kernel = {"poly0", 0, 1, &Poly0Fun};

    run_test(n, kernel, UNIF);
  }

  DEBUG("Kernel: f=xyz (Unif)");
  {
    int n = 2;
    kernel_t kernel = {"poly1", -3, -1, &Poly1Fun};

    run_test(n, kernel, UNIF);
  }

  DEBUG("Kernel: f=xyz (Cheb)");
  {
    int n = 2;
    kernel_t kernel = {"poly1", -3, -1, &Poly1Fun};

    run_test(n, kernel, CHEB);
  }

  DEBUG("Kernel: f=x^2 y^3 z (Unif)");
  {
    int n = 4;
    kernel_t kernel = {"poly3", -6,  1, &Poly3Fun};

    run_test(n, kernel, UNIF);
  }
  
  DEBUG("Kernel: f=x^2 y^3 z (Cheb)");
  {
    int n = 4;
    kernel_t kernel = {"poly3", -6,  1, &Poly3Fun};

    run_test(n, kernel, CHEB);
  }

  DEBUG("Run was successful");

  return 0;
}

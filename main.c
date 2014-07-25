#include "bbfmm.h"
#include "directCalc.h"

int main(int argc, char *argv[]) {

  /*
  // usage description
  printf("Usage:\n\t %s N L n\nwhere N is the number of source"
  " points (and field points),\n L is the length of the"
  " computational domain\n and n is the number of Chebyshev"
  " order.\n e.g.: %s 1000 1 3.\n", argv[0], argv[0]);
  */

  // default arguments of fmm
  int N          = 1000; // number of particles
  int n          = 3;    // order of the interpolation polynomial, either Chebyshev or ...
  int level_tree = 3;    // hierarchical tree level
  int level_pbc  = 0;    /* lpbc  = 0, 1, 2, 3
			  * start = 0, 1, 4, 13
			  * Note: the original computational domain is lpbc=0
			  *       and PBC starts from lpbc=2
			  */
  gridT grid     = UNIFORM;
  double box_len = 4.;   // size of simulation cell: [-len/2,len/2]^3
  double alpha   = 0.;   // adjusting box size to 1+alpha, 0<=alpha<.5. In tensor version, alpha is always 0 for best accuracy

	
  // parse input arguments
  {
    int i;
    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i],"-n"))
	N = atoi(argv[++i]);
      if (!strcmp(argv[i],"-o"))
	n = atoi(argv[++i]);
      if (!strcmp(argv[i],"-l"))
	level_tree = atoi(argv[++i]);
      if (!strcmp(argv[i],"-pbc"))
	level_pbc = atoi(argv[++i]);
      if (!strcmp(argv[i],"-g"))
	grid = atoi(argv[++i]);
      if (!strcmp(argv[i],"-s"))
	box_len = atof(argv[++i]);
      if (!strcmp(argv[i],"-a"))
	alpha = atof(argv[++i]);
    }
    printf("Starting black-box fast multipole method ... \n"
	   "Number of particles: %d,\n"
	   "Polynomial order: %d,\n"
	   "Fmm tree level: %d,\n"
	   "Fmm PBC level: %d,\n"
	   "Using grid: %d (0-unifrom; 1-Chebyshev),\n"
	   "Box size: %.2f,\n"
	   "Adjust box size by: %.2f,\n",
	   N, n, level_tree, level_pbc, grid, box_len, alpha);
  }
  


  // set problem dimension
  dof2 dof = {1, 1}; 


  // randomly set source and field points
  int Ns = N;
  int Nf = N;
  vec3 *source = malloc( Ns * sizeof(vec3) );
  vec3 *field  = malloc( Nf * sizeof(vec3) );
  double *q    = malloc( Ns * sizeof(double) );

  SetSources( field, Nf, source, Ns, q, dof.s, box_len );


  
  //AnisoInit();
  //kernelT kfun = &AnisoFun;
  //kernel_t kernel = {"lapforce", 2.0, -1, &LapforceFun};
  //kernel_t kernel = {"gaussian",  0.0,  0, &GaussianFun};
  kernel_t kernel = {"laplace", 1.0, 1, &LaplacianFun};

  
  double epsilon = 1e-9;
  //calloc( Ns, sizeof(double) );
  double *phiFmm = bbfmm( field, Nf, source, Ns, q, dof, box_len, alpha,
			  n, grid, kernel, level_tree, level_pbc,
			  epsilon );
  //print_array(phiFmm, Ns*dof.s, "FMM result");
  
        
  // Direct calculation
  // Make sure initialized to 0
  //= calloc( Nf, sizeof(double) );
  double *phiDir = directCalc( field, Nf, source, Ns, q, dof, box_len,
			       kernel.kfun, level_pbc);
  //print_array(phiDir, Nf*dof.f, "Direct calculation");
  
	
  double err = 0;
  err = ComputeL2Err(phiFmm, phiDir, Nf*dof.f);
  printf("Error of FMM: %e\n", err);

  

  //AnisoClean();
  free(field);
  free(source);
  free(q);
  free(phiFmm);
  free(phiDir);

  return 0;
}




#include "bbfmm.h"
#include "directCalc.h"

int main(int argc, char *argv[]) {

  // usage description
  printf("=======================================\n"
	 "Usage:\n"
	 "-n: number of source (field) points;\n"
	 "-o: interpolation order;\n"
	 "-l: FMM tree level;\n"
	 "-g: 0 for unifrom grids and 1 for Chebyshev grids;\n"
	 "-s: size of simulation cell ( [-length/2, length/2]^3 );\n"
	 "-pbc: FMM PBC level;\n"
	 "-a: adjusting box size to 1+alpha, 0<=alpha<.5. In tensor version, alpha is always 0 for best accuracy;\n"
	 "-dir: 1 for using direct calculation and compute accuracy and 0 for not;\n"
	 "=======================================\n\n"
	 );


  // default arguments of fmm
  int N          = 1000; // number of particles
  int n          = 3;    // order of the interpolation polynomial
  int level_tree = 3;    // hierarchical tree level
  int level_pbc  = 0;    /* lpbc  = 0, 1, 2, 3
			  * start = 0, 1, 4, 13
			  * Note: the original computational domain is lpbc=0
			  *       and PBC starts from lpbc=2
			  */
  gridT grid     = UNIFORM;
  double box_len = 4.;
  double alpha   = 0.;
  int dirCalculate = 1;
	
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
      if (!strcmp(argv[i],"-dir"))
	dirCalculate = atoi(argv[++i]);
    }
    printf("***************************************\n"
	   "Starting black-box fast multipole method ... \n"
	   "Number of particles: %d,\n"
	   "Polynomial order: %d,\n"
	   "FMM tree level: %d,\n"
	   "FMM PBC level: %d,\n"
	   "Grid: %s,\n"
	   "Box size: %.2f,\n"
	   "Adjust box size by: %.2f,\n",
	   N, n, level_tree, level_pbc, grid ? "Chebyshev grids": "uniform grids", box_len, alpha);
    printf("Direct calculation: %s.\n", dirCalculate ? "yes" : "no");
    printf("***************************************\n\n");
  }
  


  // set problem dimension
  int2 dof = {1, 1}; 


  // randomly set source and field points
  int Ns = N;
  int Nf = N;
  vec3 *source = malloc( Ns * sizeof(vec3) );
  vec3 *field  = malloc( Nf * sizeof(vec3) );
  double *q    = malloc( Ns * sizeof(double) );

  SetSources( field, Nf, source, Ns, q, dof.s, box_len );


  
  //AnisoInit();
  //kernelT kfun = &AnisoFun;

  // {name, homogen, symmetry, kernelFunction}
  kernel_t kernel = {"lapforce", 2.0, -1, &LapforceFun};
  //kernel_t kernel = {"gaussian",  0.0,  0, &GaussianFun};
  //kernel_t kernel = {"laplace", 1.0, 1, &LaplacianFun};

  
  double epsilon = 1e-9;
  double *phiFmm = bbfmm( field, Nf, source, Ns, q, dof, box_len, alpha,
			  n, grid, kernel, level_tree, level_pbc,
			  epsilon );
  //print_array(phiFmm, Ns*dof.s, "FMM result");
  
        
  // Direct calculation
  if (dirCalculate) {

    double *phiDir = directCalc( field, Nf, source, Ns, q, dof, box_len,
				 kernel.kfun, level_pbc );
    //print_array(phiDir, Nf*dof.f, "Direct calculation");  
	
    double err = ComputeL2Err(phiFmm, phiDir, Nf*dof.f);
    printf("Error of FMM: %e\n", err);
    free(phiDir);
  }

  //AnisoClean();
  free(field);
  free(source);
  free(q);
  free(phiFmm);

  return 0;
}




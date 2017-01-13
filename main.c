#include "bbfmm.h"
#include "directCalc.h"


void particleFMM(int argc, char *argv[]);
void segmentFMM (int argc, char *argv[]);

int  main       (int argc, char *argv[]) {


  //particleFMM(argc, argv);
  segmentFMM (argc, argv);
  

  return 0;
}


void particleFMM(int argc, char *argv[]) {

  // usage description
  printf("=======================================\n"
	 "Usage:\n"
	 "-n:   number of source (field) points;\n"
	 "-o:   interpolation order;\n"
	 "-l:   FMM tree level;\n"
	 "-g:   0 for unifrom grids and 1 for Chebyshev grids;\n"
	 "-s:   size of simulation cell ( [-length/2, length/2]^3 );\n"
	 "-pbc: FMM PBC level;\n"
	 "-eps: threshold to truncate singular values after the SVD;\n"
	 "-dir: 1 for using direct calculation and compute accuracy and 0 for not;\n"
	 "=======================================\n\n");


  // default arguments of fmm
  int N          = 1000; // number of particles
  int n          = 3;    // order of the interpolation polynomial
  int dirCalc    = 1;    // enable direct calculation
  int tree_lvl   = 3;    // hierarchical tree level
  int pbc_lvl    = 0;    /* lpbc  = 0, 1, 2, 3
			    * start = 0, 1, 4, 13
			    * Note: the original computational domain is lpbc=0
			    *       and PBC starts from lpbc=2
			    */  
  gridT  grid    = UNIF;
  double box_len = 4.;
  double alpha   = 0.;
  double epsilon = 1e-9;

	
  // parse input arguments
  {
    int i;
    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i],"-n"))
	N        = atoi(argv[++i]);
      if (!strcmp(argv[i],"-o"))
	n        = atoi(argv[++i]);
      if (!strcmp(argv[i],"-g"))
	grid     = atoi(argv[++i]);
      if (!strcmp(argv[i],"-pbc"))
	pbc_lvl  = atoi(argv[++i]);
      if (!strcmp(argv[i],"-l"))
	tree_lvl = atoi(argv[++i]);
      if (!strcmp(argv[i],"-dir"))
	dirCalc  = atoi(argv[++i]);
      if (!strcmp(argv[i],"-s"))
	box_len  = atof(argv[++i]);
      if (!strcmp(argv[i],"-eps"))
	epsilon  = atof(argv[++i]);
    }
    printf("***************************************\n"
	   "Starting black-box fast multipole method ... \n"
	   "Number of particles:   %d,\n"
	   "Polynomial order:      %d,\n"
	   "FMM tree level:        %d,\n"
	   "FMM PBC level:         %d,\n"
	   "Grid:                  %s,\n"
	   "Box size:              %.2f,\n",
	   N, n, tree_lvl, pbc_lvl, grid ? "Chebyshev": "uniform", box_len);
    printf("Direct calculation:    %s.\n", dirCalc ? "yes" : "no");
    printf("***************************************\n\n");
  }  


  // set problem dimension
  int2 dof = {1, 1}; 


  // randomly set source and field points
  int Ns = N;
  int Nf = N;
  vec3 *source = malloc( Ns       * sizeof(vec3) );
  vec3 *field  = malloc( Nf       * sizeof(vec3) );
  double *q    = malloc( Ns*dof.s * sizeof(double) );

  SetSources( field, Nf, source, Ns, q, dof.s, box_len );
  

  /* kernel_t struct: {name, homogen, symmetry, kernelFunction}
   * name can be arbitrary string given by user
   * homogen = log_2 ( f(r) / f(2r) )
   * symm=-1,1 for anti-symmetry and symmetry
   * kernel function has to be defined first
   */
  //kernel_t kernel = {"gaussian",  0.0,  0, &GaussianFun};
  kernel_t kernel = {"laplace", 1.0, 1, &LaplacianFun};
  //kernel_t kernel = {"poly3", 0,  1, &Poly3Fun};
  //kernel_t kernel = {"poly3", -6,  1, &Poly3Fun};
    
  
  FMMSrc fmm_src = {PNT, source, NULL, NULL, q, NULL, Ns, 0};
  double *phiFmm = bbfmm( fmm_src, field, Nf, dof, box_len, alpha,
			  n, grid, kernel, tree_lvl, pbc_lvl, epsilon );
  
        
  // Direct calculation
  if (dirCalc) {

    double *phiDir = directCalc( field, Nf, source, Ns, q, dof, box_len,
				 kernel.kfun, pbc_lvl );
	
    double err = ComputeL2Err(phiFmm, phiDir, Nf*dof.f);
    printf("Error of FMM: %e\n", err);
    free(phiDir);
  }

  free(field);
  free(source);
  free(q);
  free(phiFmm);

}


void segmentFMM(int argc, char *argv[]) {

  // usage description
  printf("=======================================\n"
	 "Usage:\n"
	 "-o:   interpolation order;\n"
	 "-l:   FMM tree level;\n"
	 "-g:   0 for unifrom grids and 1 for Chebyshev grids;\n"
	 "-s:   size of simulation cell ( [-length/2, length/2]^3 );\n"
	 "-a:   adjusting box size to 1+alpha, 0<=alpha<.5. In tensor version, alpha is always 0 for best accuracy;\n"
	 "-pbc: FMM PBC level;\n"
	 "-eps: threshold to truncate singular values after the SVD;\n"
	 "-dir: 1 for using direct calculation and compute accuracy and 0 for not;\n"
	 "=======================================\n\n");


  // default arguments of fmm
  int Ns         = 3;
  int Nf         = 1;
  int n          = 3;    // order of the interpolation polynomial
  int dirCalc    = 1;    // enable direct calculation
  int tree_lvl   = 3;    // hierarchical tree level
  int pbc_lvl    = 0;    /* lpbc  = 0, 1, 2, 3
			  * start = 0, 1, 4, 13
			  * Note: the original computational domain is lpbc=0
			  *       and PBC starts from lpbc=2
			  */
  gridT  grid    = UNIF;
  double box_len = 4.;
  double alpha   = 0.;
  double epsilon = 1e-9;
	
  // parse input arguments
  {
    int i;
    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i],"-o")) {
	n        = atoi(argv[++i]);
	epsilon  = pow(10,-n); // set Chebyshev compression ratio
      }
      if (!strcmp(argv[i],"-l"))
	tree_lvl = atoi(argv[++i]);
      if (!strcmp(argv[i],"-pbc"))
	pbc_lvl  = atoi(argv[++i]);
      if (!strcmp(argv[i],"-g"))
	grid     = atoi(argv[++i]);
      if (!strcmp(argv[i],"-dir"))
	dirCalc  = atoi(argv[++i]);
      if (!strcmp(argv[i],"-s"))
	box_len  = atof(argv[++i]);
      if (!strcmp(argv[i],"-a"))
	alpha    = atof(argv[++i]);
      if (!strcmp(argv[i],"-eps"))
	epsilon  = atof(argv[++i]);
    }
    printf("***************************************\n"
	   "Starting black-box fast multipole method ... \n"
	   "Number of segments:    %d,\n"
	   "Number of field point: %d,\n"
	   "Polynomial order:      %d,\n"
	   "FMM tree level:        %d,\n"
	   "FMM PBC level:         %d,\n"
	   "Grid:                  %s,\n"
	   "Box size:              %.2f,\n"
	   "Adjust box size by:    %.2f,\n",
	   Ns, Nf, n, tree_lvl, pbc_lvl, grid ? "Chebyshev": "uniform", box_len, alpha);
    printf("Direct calculation:    %s.\n", dirCalc ? "yes" : "no");
    printf("***************************************\n\n");
  }

  // set problem dimension {source, field}
  int2 dof = {9, 6}; 


  // Segments are the three edges of a triangle
  // note: be careful of the box size, use box_len=4 (from -1 to 1) in this case 
  segT segment[Ns];
  segment[0] = (segT){ (vec3){0.620000, 0.590000, 0.570000},
		       (vec3){0.710000, 0.510000, 0.650000} };
  segment[1] = (segT){ (vec3){0.710000, 0.510000, 0.650000},
		       (vec3){0.900000, 0.510000, 0.650000} };
  segment[2] = (segT){ (vec3){0.900000, 0.510000, 0.650000},
		       (vec3){0.620000, 0.590000, 0.570000} };

 
  // Burger's vector
  double burger[] = { -1.190000, -0.780000, -0.330000,
		      -1.190000, -0.780000, -0.330000,
		      -1.190000, -0.780000, -0.330000, };
  
  // one field point
  vec3 field[] = { (vec3){-1.42, 1.37, 1.72} };
  
  
  /* kernel_t struct: {name, homogen, symmetry, kernelFunction}
   * name can be arbitrary string given by user
   * homogen = log_2 ( f(r) / f(2r) )
   * kernel function has to be defined first
   */
  AnisoInit();
  //kernel_t kernel = {"lapforce", 2.0, -1, &LapforceFun};
  kernel_t kernel = {"Aniso", 2.0, -1, &AnisoFun};
  
  FMMSrc fmm_src = {SEG, NULL, NULL, segment, NULL, burger, Ns, n+2};
  double *phiFmm = bbfmm( fmm_src, field, Nf, dof, box_len, alpha,
  			  n, grid, kernel, tree_lvl, pbc_lvl, epsilon );

#if 0
  printf("FMM with lpbc=%d:\n", pbc_lvl);
  int i;
  for (i=0; i<Nf*dof.f; i++)
    printf("%e\t", phiFmm[i]);
  printf("\n");
#endif
    
  // Direct calculation
  if (dirCalc) {

    int    srcSize = fmm_src.Ns * fmm_src.nGauss;
    vec3   seg2pnt[srcSize];
    double q[srcSize * dof.s];
    Segment2Point(fmm_src.segment,
		  fmm_src.Ns,
		  fmm_src.burger,
		  fmm_src.nGauss,
		  seg2pnt,
		  q);

    // Check zero net density
    int i;
    double sum = 0;
    for (i=0; i<srcSize*dof.s; i++)
      sum += q[i];
    assert(sum < 1e-12);

      
    // Direct calculation
    double *phiDir = directCalc( field, Nf, seg2pnt, srcSize, q, dof, box_len,
				 kernel.kfun, pbc_lvl );

#if 0
    printf("Direct with lpbc=%d:\n", pbc_lvl);
    for (i=0; i<Nf*dof.f; i++)
      printf("%e\t", phiDir[i]);
    printf("\n");
#endif

    double err = 0;
    err = ComputeL2Err(phiFmm, phiDir, Nf*dof.f);
    printf("Error of FMM with lpbc=%d: %e\n", pbc_lvl, err);
 
    free(phiDir), phiDir=NULL;
  }

  AnisoClean();
  free(phiFmm);

}


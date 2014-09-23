#include <assert.h>   // assert()
#include <fftw3.h>    // fft transform of real input
#include <math.h>

#include "bbfmm.h"


#define FFTW_FLAG       FFTW_ESTIMATE // option for fftw plan type
#define HOMO_THRESHOLD  1e-1          // threshold for homogeneous kenrel

//#define NUMGAUSS 10                   // Gauss quadrature points
// choose to match the interpolation accuracy

  
double* bbfmm( FMMSrc fmm_src, vec3 *field, int Nf, int2 dof,
	       double boxLen, double alpha, int n, gridT grid_type,
	       kernel_t kernel, int treeLevel, int lpbc,
	       double epsilon ) {

  
  vec3 *source = fmm_src.source;
  int   Ns     = fmm_src.Ns;

    
  timeType t0 = Timer();          // Begin pre-computation timing


  /* ---- Compute the Chebyshev weights and the lookup table ---- */

  int n2 = n*n;          // n2 = n^2
  int n3 = n2*n;         // n3 = n^3

  int    Ktable  [ 343]; // Cell index within three layers on every side: 7^3 = 343
  double Kweights[  n3]; // Omega matrix for SVD
  double Cweights[2*n2]; // X2X operator (Sn on page 8715)
  double Tkz     [  n2]; // Evaluation of n chebyshev nodes (of T_n) on
			 // chebyshev polynomials from T_0 to T_{n-1}

  ComputeWeights(Tkz,Ktable,Kweights,Cweights,n,alpha,grid_type);


  /* ---- Build FMM tree ---- */  

  nodeT  *tree = BuildTree(boxLen, treeLevel, n);


  /* ---- Precompute M2L operator ---- */
    
  int2    cutoff;     // Compression index of SVD
  double *K  = NULL;  // K: C^{(i)} the compressed M2L operator 
  double *U  = NULL;  // U: U^k_r p. 8719; downward pass; field
  double *VT = NULL;  // V: S^K_r p. 8718; upward pass; source

  GetM2L(Kweights, boxLen, alpha, &cutoff, n, kernel, epsilon, dof,
	 treeLevel, &K, &U, &VT, grid_type);
    

  if (fmm_src.src_t == SEG) {
    fmm_src.midpnt  = ComputeMiddlePoint( fmm_src.segment, Ns );
    FMMDistribute(tree, field, fmm_src.midpnt, Nf, Ns, treeLevel);
  } else 
    FMMDistribute(tree, field, source,         Nf, Ns, treeLevel);

  
  // End pre-computation timing
  timeType t1 = Timer();


  /* ---- FMM computation ---- */
  
  // allocate the FMM result array
  double *stressFmm = FMMCompute(&tree, fmm_src, field, Nf, K, U, VT, Tkz,
				 Ktable, Kweights, Cweights, cutoff, n, dof,
				 alpha, grid_type, boxLen, lpbc, kernel);
  
	
  // End FMM timing
  timeType t2 = Timer();
    	

  /* ---- Free the resources ---- */
  
  FMMCleanup(tree);

  free(K);
  if (U  != NULL)
    free(U );
  if (VT != NULL)
    free(VT);

  if (fmm_src.midpnt != NULL)
    free(fmm_src.midpnt);


  /* ---- Output and return ---- */
  
  printf("Pre-Comp Time : %f seconds\n", t1-t0);
  printf("FMM      Time : %f seconds\n", t2-t1);

  return stressFmm;
}


// Builds the FMM hierarchy
nodeT* BuildTree(double boxLen, int treeLevel, int n) {

  vec3   center = {0, 0, 0};
  nodeT* root   = NewNode(center, boxLen, n);
  BuildFMMHierarchy(root, treeLevel, n);

  return root;
}

 
/*
 * Function: ComputeM2L
 * -----------------------------------------------------------------
 * Prepare for the FMM calculation by pre-computing the SVD (if necessary),
 * and reading in the necessary matrices.
 */
void GetM2L(double *Kweights, double boxLen, double alpha,
	    int2 *cutoff, int n, kernel_t kernel,
	    double epsilon, int2 dof, int treeLevel,
	    double **K, double **U, double **VT, int grid_type) {

   
  char *dir_name = "precompute/"; 
  create_directory( dir_name );

  char Kmat[50];
  char Umat[50];
  char Vmat[50];
  char *grid1 = "Cheb";
  char *grid2 = "Unif";

   
  if (grid_type == CHEB) {
    sprintf(Kmat,"%s%s%sK%da%.1f.bin", dir_name, kernel.name, grid1, n, alpha);
    sprintf(Umat,"%s%s%sU%da%.1f.bin", dir_name, kernel.name, grid1, n, alpha);
    sprintf(Vmat,"%s%s%sV%da%.1f.bin", dir_name, kernel.name, grid1, n, alpha);
  } else {
    sprintf(Kmat,"%s%s%sK%da%.1f.bin", dir_name, kernel.name, grid2, n, alpha);
    sprintf(Umat,"Readme.txt"); // uniform grid does not have U or V file,
    sprintf(Vmat,"Makefile");   // so just make sure thest two files exist
  }
	

  if ( !PrecomputeAvailable(Kmat, Umat, Vmat, kernel.homogen, boxLen, treeLevel, grid_type) ) {
    StartPrecompute( boxLen, treeLevel, n, dof, kernel, Kmat, Umat, Vmat, alpha, Kweights, epsilon, grid_type );
  }
   
  // Read kernel interaction matrix K and singular vectors U and VT
  FMMReadMatrices(K, U, VT, cutoff, n, dof, Kmat, Umat, Vmat,
		  treeLevel, kernel.homogen, grid_type);

}


// check if the precompute files exist or if the existing files are usable
bool PrecomputeAvailable( char *Kmat, char *Umat, char *Vmat,
			  double homogen, double boxLen,
			  int treeLevel, int grid_type ) {
  
  FILE *fK, *fU, *fV;
  fK = fopen(Kmat, "rb");
  fU = fopen(Umat, "rb");
  fV = fopen(Vmat, "rb");
	 
  if (fK == NULL || fU == NULL || fV == NULL) { // files do not exist
    
    if (fK!=NULL) fclose(fK);
    if (fU!=NULL) fclose(fU);
    if (fV!=NULL) fclose(fV);
    return false;
    
  } else if ( !IsHomoKernel(homogen) ) { // non-homogeneous kernel

    double len_file;
    int    lvl_file;    
    READ_CHECK( fread(&len_file, sizeof(double), 1, fK), 1 );
    READ_CHECK( fread(&lvl_file, sizeof(int),    1, fK), 1 );
    fclose(fK);
    fclose(fU);
    fclose(fV);
    
    if ( lvl_file != treeLevel || !equal_double( len_file, boxLen ) )   
      return false;
    else
      return true;    
  }

  fclose(fK);
  fclose(fU);
  fclose(fV);
  return true;
}


void StartPrecompute(double boxLen, int treeLevel, int n, int2 dof, kernel_t kernel, char *Kmat, char *Umat, char *Vmat, double alpha, double *Kweights, double epsilon, int grid_type) {

  printf("Generate precompute file ...\n");
  FILE *fK = fopen(Kmat, "wb");
  fwrite(&boxLen, sizeof(double), 1, fK); // write box size
  fwrite(&treeLevel, sizeof(int), 1, fK); // write tree level
  fclose(fK);	     		 

  if ( IsHomoKernel(kernel.homogen) ) {  // homogeneours kernel

    double unit_len = 1.0;
    compute_m2l_operator(n, dof, kernel, Kmat, Umat, Vmat, unit_len, alpha, Kweights, epsilon, grid_type);
    
  } else {                                // non-homogeneous kernel

    int j;
    double boxLenLevel = boxLen/4;        // FMM starts from the second level
    for (j=2; j<=treeLevel; j++) {
      compute_m2l_operator(n, dof, kernel, Kmat, Umat, Vmat, boxLenLevel, alpha, Kweights, epsilon, grid_type);
      boxLenLevel/=2;
    }

  }                                       // end non-homogeneous kernel
}


void compute_m2l_operator (int n, int2 dof, kernel_t kernel, char *Kmat, char *Umat, char *Vmat, double l, double alpha, double *Kweights, double epsilon, int grid_type) {

  switch (grid_type) {
      
  case UNIF:
    ComputeKernelUnif(n, dof, kernel.kfun, Kmat, alpha, l);
    break;
  case CHEB:
    ComputeKernelCheb(Kweights, n, kernel, epsilon, dof,
		      Kmat, Umat, Vmat, alpha, l);
    break;
  }
}

 
/*
 * Function: FMMReadMatrices
 * ------------------------------------------------------------------
 * Read in the kernel interaction matrix M and the matrix of singular
 * vectors U.
 */
void FMMReadMatrices(double **K, double **U, double **VT, int2 *cutoff,
		     int n, int2 dof, char *Kmat, char *Umat,
		     char *Vmat, int treeLevel, double homogen,
		     int grid_type) {

  int preCompLevel;
  if ( !IsHomoKernel(homogen) )
    preCompLevel = treeLevel - 1;
  else
    preCompLevel = 1;
  

  unsigned long long Ksize;
  int Usize;
  int Vsize;
  
  if (grid_type == UNIF) {
       
    cutoff->f = dof.f*n*n*n;
    cutoff->s = dof.s*n*n*n;

    Ksize = 316*(2*n-1)*(2*n-1)*(2*n-1)*dof.s*dof.f * preCompLevel;
    
  } else if (grid_type == CHEB) { // read 'U' and 'V'
    
    FILE *fU = fopen(Umat,"rb"); 
    FILE *fV = fopen(Vmat,"rb");

    READ_CHECK( fread(&(cutoff->f), sizeof(int), 1, fU), 1 );
    READ_CHECK( fread(&(cutoff->s), sizeof(int), 1, fV), 1 );

    Usize = cutoff->f * dof.f * n*n*n * preCompLevel;
    Vsize = cutoff->s * dof.s * n*n*n * preCompLevel;
    Ksize = 316 * cutoff->f*cutoff->s * preCompLevel;
    
    *U  = (double *)malloc(Usize *sizeof(double)); 
    *VT = (double *)malloc(Vsize *sizeof(double));
    
    READ_CHECK( fread(*U,  sizeof(double), Usize, fU), Usize );
    READ_CHECK( fread(*VT, sizeof(double), Vsize, fV), Vsize );

    fclose(fU);
    fclose(fV);
    
  } else
    assert(false);

  
  FILE *fK  = fopen(Kmat, "rb");
  assert(Ksize > 0); // check over flow
  *K = malloc(Ksize *sizeof(double));    
  fseek(fK, 1*sizeof(int) + 1*sizeof(double), SEEK_SET); // skip 'tree level' and 'box length'
  READ_CHECK( fread(*K, sizeof(double), Ksize, fK), Ksize );
  fclose(fK);

}

 
/*
 * Function: FMMDistribute
 * ------------------------------------------------------------------
 * Distribute the field points and sources to the appropriate location
 * in the FMM hierarchy and sets up the interaction list.
 */
void FMMDistribute(nodeT *A, vec3 *field, vec3 *source, int Nf, 
		   int Ns, int level) {

  int *fieldlist  = (int *)malloc(Nf * sizeof(int));
  int *sourcelist = (int *)malloc(Ns * sizeof(int));  
	
  // Initialize the point distribution for the root node
  int i;
  for (i=0;i<Nf;i++)
    fieldlist[i]  = i;
  for (i=0;i<Ns;i++)
    sourcelist[i] = i;
  A->Nf = Nf;
  A->Ns = Ns;
     
  // Distribute all of the source and field points to the appropriate cell
  DistributeFieldPoints(A,field,fieldlist,level);
  DistributeSources(A,source,sourcelist,level);
	
  // Construct the interaction list for the entire octree
  A->neighbors[0] = A;
  A->ineigh = 1;
  A->cshiftneigh[0] = (vec3){0, 0, 0};
  InteractionList(A,level);  // build the interaction that need to be computed
	
  free(fieldlist);
  free(sourcelist);
}

/*
 * Function: FMMCompute
 * ------------------------------------------------------------------
 * Computes the field using BBFMM.
 */
/*
  #ifdef LINEINT
  void FMMCompute(nodeT **A, FMMSrc fmm_src, vec3 *field, int Nf, segT *segment, vec3 *
  midpoint, double *burg, double *K, double *U,
  double *VT, double *Tkz, int *Ktable,
  double *Kweights, double *Cweights,
  int2 cutoff, int n, int2 dof, double alpha,
  double *phi, int grid_type, double boxLen, int lpbc,
  kernel_t kernel) {
    
  // Compute the Guassian quadrature points and weights
  double *gpoint  = (double *) malloc(NUMGAUSS * sizeof(double));
  double *gweight = (double *) malloc(NUMGAUSS * sizeof(double));
  gqwp(numgau, gpoint, gweight);
    

  printf("Begin upward pass ...\n");
  UpwardPass(A, segment, midpoint, Cweights, Tkz, burg, n,
  dof.s, gpoint, gweight, alpha, grid_type);
    
    
  free(gpoint);
  free(gweight);
    
  #elif TENSOR
*/
     
double* FMMCompute(nodeT **A, FMMSrc fmm_src, vec3 *field, int Nf,
		   double *K, double *U, double *VT, double *Tkz,
		   int *Ktable, double *Kweights, double *Cweights,
		   int2 cutoff, int n, int2 dof, double alpha, 
		   int grid_type, double boxLen, int lpbc,
		   kernel_t kernel) {

  
  if (fmm_src.src_t == SEG)
    InitGaussQuadrature(fmm_src.nGauss);

  printf("Begin upward pass ...\n");
  timeType t1 = Timer();
  UpwardPass(A, fmm_src, Cweights, Tkz, n, dof.s, alpha, grid_type);
  timeType t2 = Timer();
  printf("Time for upward pass: %f\n", t2-t1);
  
  if (fmm_src.src_t == SEG)
    CleanGaussQuadrature();


  double *phi = calloc( Nf*dof.f, sizeof(double) );
  AddPBC( lpbc, (*A)->sourceval, &((*A)->fieldval), phi, Nf,
	  n, dof, boxLen, alpha, Tkz, kernel, grid_type );

  // Compute cell interactions.
  printf("Begin M2L ...\n");
  int root_level = 0;
  t1 = Timer();
  FMMInteraction(A, K, Ktable, U, VT, Kweights, n, dof, cutoff,
		 kernel.homogen, root_level, grid_type);  
  t2 = Timer();
  printf("Time for upward pass: %f\n", t2-t1);

  printf("Begin downward pass ...\n");
  t1 = Timer();
  DownwardPass(A, field, fmm_src, Cweights, Tkz, n, dof,
	       alpha, kernel.kfun, phi, grid_type);
  t2 = Timer();
  printf("Time for upward pass: %f\n", t2-t1);

  
  return phi; 
}


/*
 * Function: FMMCleanup
 * ------------------------------------------------------------------
 * Cleans up after FMM calculation.
 */
void FMMCleanup(nodeT *A) {

  FreeNode(A);
  fftw_cleanup();
}

/*
 * Function: ReadParam
 * ------------------------------------------------------------------
 * Read in the user-specified options from the options file.
 * 
 */
// TODO: read cutoff and dof
/*
  int SetParam(char *options, fmmparam *param) {
  char tag[10], line[80];
  int i, tag_id, ntags=13;
  char alltags[][10]={"Ns","Nf","dof","L","n","levels","PBClevels",
  "PBCshells","precomp","cutoff","homogen",
  "filesval","filesvec"};
  enum {Ns,Nf,dof,L,n,levels,PBClevels,PBCshells,precomp,cutoff,homogen,
  filesval,filesvec};
	
  FILE *f;
	
  // Initializes the parameter values
  param->Ns = 10;
  param->Nf = 10;
  param->dof = 3;
  param->L    = 1;
  param->n    = 4;
  param->levels = 3;
  param->PBClevels  = 0;
  param->PBCshells = 0;
  param->precomp = 1;
  param->cutoff = 96;
  param->homogen = 1;
  //param->filesval="t.out";
  //param->filesvec="u.out";
	
  // Open options file for reading
  f = fopen(options,"r");
	
  // Reads one line at a time until EOF is reached
  while (fgets(line,80,f) != NULL) {
  // Read in tag
  sscanf(line,"%s",tag);
		
  // Determine the id corresponding to the tag
  tag_id = -1;
  for (i=0;i<ntags;i++) {
  if (strcmp(tag,alltags[i]) == 0) {
  tag_id = i;
  break;
  }
  }
		
  // Stores the selected parameter in the appropriate variable
  switch (tag_id) {
  // Read in number of sources
  case Ns:
  sscanf(line,"%s %d",tag,&param->Ns);
  break;
				
  // Read in number of field points
  case Nf:
  sscanf(line,"%s %d",tag,&param->Nf);
  break;
				
  // Read in number of degrees of freedom : 
  case dof:
  sscanf(line,"%s %d",tag,&param->dof);
  break;
				
  // Read in length of one side of simulation cube
  case L:
  sscanf(line,"%s %lf",tag,&param->L);
  break;
				
  // Read in number of Chebyshev nodes in one direction
  case n:
  sscanf(line,"%s %d",tag,&param->n);
  break;
				
  // Read in maximum number of levels in octree
  case levels:
  sscanf(line,"%s %d",tag,&param->levels);
  break;
				
  // Read in number of PBC levels
  case PBClevels:
  sscanf(line,"%s %d",tag,&param->PBClevels);
  break;
				
  // Read in number of PBC shells
  case PBCshells:
  sscanf(line,"%s %d",tag,&param->PBCshells);
  break;
				
  // Read in pre-computation switch  -> not used
  case precomp:
  sscanf(line,"%s %d",tag,&param->precomp);
  break;
				
  // Read in number of singular values to keep : max value is n^3. Try to reduce it to adjust the error
  case cutoff:
  sscanf(line,"%s %d",tag,&param->cutoff);
  break;
				
  // Read in order of homogeneity of kernel  1/R = 1 1/R^2 -> 2
  case homogen:
  sscanf(line,"%s %lf",tag,&param->homogen);
  break;
				
  // Read in file name for singular values -> not used
  case filesval:
  sscanf(line,"%s %s",tag,param->filesval);
  break;
				
  // Read in file name for singular vectors -> not used
  case filesvec:
  sscanf(line,"%s %s",tag,param->filesvec);
  break;
				
  default:
  printf("%s is an invalid tag!\n",tag);
  }
  }
	
  // Closes options file
  fclose(f);
	
  return 0;
  }
*/


/* Input:
 * 	sourcelist - index of the source (segment) point
 * 	num - number of gaussian points used
 * Output:
 * 	sourcet - Gauss points along the segments transformed into a unit cube
 *	xi - Jacobi of the integral
 */ 
void LineIntegral( FMMSrc fmm_src, int *sourcelist, int Ns, double L,
		   vec3 scenter, vec3 *sourcet, double *xi){
  

  int   nGauss   = fmm_src.nGauss;
  segT *segment  = fmm_src.segment;
  vec3 *midpoint = fmm_src.midpnt;
  
  vec3 *source = (vec3 *) malloc(nGauss * Ns * sizeof(vec3));
  vec3 p1, p2, midpt;

  double lx, ly, lz;
  int i, j, k, count;

  count = 0;
  for (i=0; i<Ns; i++) {
    k  = sourcelist[i];
    p1 = segment[k].p_beg;
    p2 = segment[k].p_end;
    midpt = midpoint[k];

    // Check: the sign cosistant with the order of gauss points
    xi[3*i + 0] = (p2.x - p1.x) / 2;
    xi[3*i + 1] = (p2.y - p1.y) / 2;
    xi[3*i + 2] = (p2.z - p1.z) / 2;

    lx = xi[3*i+0];
    ly = xi[3*i+1];
    lz = xi[3*i+2];

    // Gaussian points along the segment
    for (j=0; j<nGauss; j++, count++) {
      source[count].x = midpt.x + GAUSSP[j] * lx;
      source[count].y = midpt.y + GAUSSP[j] * ly;
      source[count].z = midpt.z + GAUSSP[j] * lz;
    }
  }

  // Map all of the source points to the box ([-1 1])^3
  double ihalfL = 2.0 / L;
  int numSource = Ns * nGauss;

  for (j=0;j<numSource;j++) {
    sourcet[j].x = ihalfL*(source[j].x - scenter.x);
    sourcet[j].y = ihalfL*(source[j].y - scenter.y);
    sourcet[j].z = ihalfL*(source[j].z - scenter.z);

    /*
      InBoxCheck(sourcet[j].x);
      InBoxCheck(sourcet[j].y);
      InBoxCheck(sourcet[j].z);
    */
  }

  free(source);		
}
   

/*
 * Function: EvaluateKernelCell
 * -------------------------------------------------------------------
 * Evaluate the kernel for interactions between a pair of cells.
 */
void EvaluateKernelCell(vec3 *field, vec3 *source, int Nf, 
			int Ns, int2 dof, kfun_t kfun,
			double *kernel) {
  
  int i, j, k, l, count, count_kernel;
  int dof_f = dof.f;
  int dof_s = dof.s;
  int LDA_kernel = dof.f * Nf;
  int dof2  = dof.f * dof.s;
  double Kij[dof2];
      
	
  for (j=0;j<Ns;j++) {
    for (i=0;i<Nf;i++) {

      kfun(&field[i], &source[j], Kij);
			
      count_kernel = dof_f * i + LDA_kernel * dof_s * j;
      count = 0;			
      for (k=0;k<dof_s;k++)
	for (l=0;l<dof_f;l++, count++)
	  // Column-major storage
	  kernel[count_kernel + k * LDA_kernel + l] = Kij[count];
    }
  }

}

/*
 * Function: EvaluateField
 * -------------------------------------------------------------------
 * Evaluate the field due to interactions between a pair of cells.
 * P2P kernel.
 * Compute K*q
 * K is (dof->f*Nf) x (dof->s*Ns)
 * q is a vector of length dof->s*Ns
 * q stores q[Ns][dof->s] in row major
 * so the result fieldval[Nf][dof->f] is also in row major
 */
void EvaluateField(vec3 *field, vec3 *source, double *q, int Nf,
		   int Ns, int2 dof, kfun_t kfun, double *fieldval) {
          

  int i, j, k, l, count, count_kernel;
  int dof2 = dof.f * dof.s;
  int LDA  = dof.f * Nf;
  int N    = dof.s * Ns;
   
  double Kij[dof2];                                 // the point to point kernel
  double* Kcell = malloc(LDA * N * sizeof(double)); // Kcell stores all the Kij's

  for (j=0;j<Ns;j++) {
    for (i=0;i<Nf;i++) {

      kfun(&field[i], &source[j], Kij);

       
      // check kernel singularity
      if ( equal_vec3(&field[i], &source[j]) )	 
	for (k=0;k<dof2;k++)
	  if ( isinf(Kij[k]) )
	    Kij[k] = 0;
       

      count_kernel = dof.f * i + LDA * dof.s * j;
      count = 0;			
      for (k=0;k<dof.s;k++)
	for (l=0;l<dof.f;l++, count++) 
	  // Kcell and Kij are column-major storage 
	  Kcell[count_kernel + k * LDA + l] = Kij[count]; 
    }
  }
	
  // Do matrix vector product
  double alpha = 1, beta = 0;
  int incr = 1;
  char trans = 'n';

  dgemv_(&trans, &LDA, &N, &alpha, Kcell, &LDA, q, &incr, &beta, fieldval, &incr);
    
  free(Kcell);
}

/*
 * Function: ComputeWeights
 * ------------------------------------------------------------------
 * Compute the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) and set up the lookup table.
 */
void ComputeWeights(double *Tkz, int *Ktable, double *Kweights, 
		    double *Cweights, int n, double alpha, int grid_type) {

  int  n3 = n*n*n;    // n3 = n^3
  int  Nc = 2*n3;     // Number of child Chebyshev nodes
  vec3 fieldt[Nc];   // Chebyshev-transformed coordinates
  vec3 Sn[n*Nc];

  // lookup table for M2L
  Create_lookup_table(Ktable);
      
  ComputeTkz( Tkz, n );
      
  // Weights for the SVD compressed matrix
  int    l1, l2, l3;
  int    count = 0;
  double tmp1;
  double tmp2;
  double nodes[n];
  GridPos1d(0, 1.0, n, grid_type, nodes);
  if (grid_type == CHEB) {
    for (l1=0;l1<n;l1++) {
      tmp1 = 1./sqrt(1-nodes[l1]*nodes[l1]);
      for (l2=0;l2<n;l2++) {
	tmp2 = tmp1/sqrt(1-nodes[l2]*nodes[l2]);
	for (l3=0;l3<n;l3++) {
	  Kweights[count] = tmp2/sqrt(1-nodes[l3]*nodes[l3]);
	  count++;
	}
      }
    }
  }

  /*
    else if (grid_type == UNIF) { // Uniform weights

    for (count=0; count<n3; count++)
    Kweights[count] = 1;
    }
  */
     
  // Map all Chebyshev nodes from the children cells to the parent
  // Note: 'nodes' below corresponds to alpha=0
  // scaled child grids: nodes[i]*(1+alpha) +- 1
  // scaled parent box: (1+alpha)*4
  // ihalfL: 1/((1+alpha)*2)
  // so the relative position is: (nodes[i] +- 1/(1+alpha)) / 2
  int    i;
  int    k    = 0;
  double L    = AdjustBoxSize(1, alpha);
  double half = 1.0/2.0; // the parent box is divided into two children in each dimension
  vec3   vtmp = {-1, -1, 0};
  for (i=0;i<2;i++) {

    // mapping function for each child cell
    if (i == 0)
      vtmp.z = -1;
    else
      vtmp.z = 1;
		
    for (l1=0;l1<n;l1++) {
      for (l2=0;l2<n;l2++) {
	for (l3=0;l3<n;l3++) {
	  fieldt[k].x = (nodes[l1] + vtmp.x / L) * half;
	  fieldt[k].y = (nodes[l2] + vtmp.y / L) * half;
	  fieldt[k].z = (nodes[l3] + vtmp.z / L) * half;
	  k++;
	}
      }
    }
  }
      
  // M2M and L2L operators
  ComputeSn(fieldt, Nc, Sn, n, Tkz, grid_type);
	
  // Extract out the Chebyshev weights
  count = 0;
  for (l1=0;l1<n;l1++) {
    k = l1*Nc;
    for (l2=0;l2<n;l2++) {
      Cweights[count] = Sn[k+l2].z;
      count++;
    }
  }
  for (l1=0;l1<n;l1++) {
    k = l1*Nc;
    for (l2=0;l2<n;l2++) {
      Cweights[count] = Sn[k+n3+l2].z;
      count++;
    }
  }
}

     
/*
 * Function: ComputeKernelCheb
 * -------------------------------------------------------------
 * Computes the kernel for 316n^6 interactions between Chebyshev nodes
 * and then computes the SVD of the kernel matrix.
 * Note: this routine utilizes the anti-symmetry of the kernel
 */
void ComputeKernelCheb(double *Kweights, int n, kernel_t
		       kernel, double epsilon, int2 dof, char
		       *Kmat, char *Umat, char *Vmat,
		       double alphaAdjust, double boxLen) {
   
  kfun_t kfun = kernel.kfun;
  int symm    = kernel.symm;
      
  static int callTime = -1;
  callTime++; // callTime = 0 for the first time called

  int i, j, l, m, m1, k1, k2, k3, l1, z;
  int count, count1, count2, count3;
  double sweight, fweight;
	
  int n3 = n*n*n;            // n3 = n^3
  int dofn3_s = dof.s * n3;
  int dofn3_f = dof.f * n3;
  int dof2n6 = dofn3_s * dofn3_f; // Total size
  int Sigma_size;
  int2 cutoff;
	
  double *K0, *U0, *Sigma, *VT0;
  double *ker, *work;
  vec3 fieldpos[n3], sourcepos[n3];
	
  K0  = (double *)malloc(316 * dof2n6 * sizeof(double));
  ker = (double *)malloc(dof2n6 * sizeof(double));

  // 316 M2L operators
  double* Kcell[316];
  for (z = 0; z < 316; ++z)
    Kcell[z] = (double *) malloc(dof2n6 * sizeof(double));

      
  // Compute the locations of the field points in a unit cube
  vec3 center = {0, 0, 0};
  GridPos3d(center, alphaAdjust, boxLen, n, CHEB, fieldpos);
  //GetBoxNodes(fieldpos, n, boxLen, center, CHEB);
      
  // Compute the kernel values for interactions with all 316 cells
  int countM2L=0;
      
  int symmNum = 158*(2-abs(symm)); // symmNum = 158, 316 for symm=pm1, 0
  int col, row, idx1, idx2;
  count=0;
  for (k1=-3;k1<4;k1++) {
    center.x = (double)k1;
    for (k2=-3;k2<4;k2++) {
      center.y = (double)k2;
      for (k3=-3;k3<4;k3++) {
	center.z = (double)k3;

	vec3 r = {k1, k2, k3};
	if ( Is_well_separated(r, 1.0) ) {
	  if (countM2L < symmNum) {

	    // Compute the source points in the cell
	    GridPos3d(center, alphaAdjust, boxLen, n, CHEB, sourcepos);
	    //GetBoxNodes(sourcepos, n, boxLen, center, CHEB);
              
	    // Compute the kernel at each of the field points 
	    EvaluateKernelCell(fieldpos, sourcepos, n3, n3, dof,
			       kfun, ker);
			 
	    // Copy the kernel values to the appropriate location
	    count1=0, count2=0;
	    for (l1=0; l1<n3; l1++, count2++) {
	      sweight = Kweights[count2];
	      for (l=0;l<dof.s;l++) { 
		for (count3=0, m1=0; m1<n3; m1++, count3++) {
		  fweight = Kweights[count3];
		  for (m=0; m<dof.f; m++, count++, count1++) { 
		    K0[count] = ker[count1]/(sweight*fweight);
		  }
		}
	      }
	    }
	    countM2L++;
	  }
	  else { // Use kernel symmetry
	    for (col=0; col<dofn3_s; col++)
	      for (row=0; row<dofn3_f; row++) {
                  
		idx1 = (315-countM2L)*dof2n6 + col*dofn3_f + row;
		idx2 = countM2L*dof2n6 +
		  (row/dof.f*dof.s + col%dof.s) *dofn3_f
		  + (col/dof.s*dof.f + row%dof.f);

		// symm=-1,1 for anti-symmetry and symmetry
		K0[idx2] = symm * K0[idx1];
		//K0[idx2] = K0[idx1];
	      }
	    countM2L++;
	  }
	}
      }
    }
  }
    
	
  // Extract the submatrix for each of the 316 cells
  count = 0;
  for (i=0;i<316;i++) {
    for (j=0;j<dof2n6;j++) {
      Kcell[i][j] = K0[count];
      count++;
    }
  }
	

  /****
   * Compute the SVD of K_fat
   ****/
	
  // Compute the SVD of K0
  char save[]="S", nosave[]="N";
  int nosavedim=1;
  int info, lwork;
  int cols_s = 316*dofn3_s;
	
  /* See dgesvd documentation:
   *          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
   *             - PATH 1  (M much larger than N, JOBU='N')  - our case for K_thin
   *             - PATH 1t (N much larger than M, JOBVT='N') - our case for K_fat
   */
  int max_dof = dof.s > dof.f ? dof.s : dof.f;
  lwork = 5*max_dof*n3; // Change
  work  = (double *) malloc(lwork * sizeof(double));
	
  int U0size;
  Sigma_size = dofn3_f;	
  U0size     = dofn3_f * dofn3_f;
  Sigma      = (double *)malloc(Sigma_size * sizeof(double));
  U0         = (double *)malloc(U0size * sizeof(double));	
  VT0        = NULL;

  dgesvd_(save, nosave, &dofn3_f, &cols_s, K0, &dofn3_f, Sigma, U0,
	  &dofn3_f, VT0, &nosavedim, work, &lwork, &info);

  FILE *ptr_file;

  double sumsqr, sum, epsqr;
  if (callTime == 0) { // The first time called
	 
    // Determine the number of singular values to keep. We use epsilon for this.
    sumsqr = 0;
    for (i=Sigma_size-1; i>=0; --i)
      sumsqr += Sigma[i] * Sigma[i];
	 
    sum = 0, epsqr = sumsqr * epsilon * epsilon;
    cutoff.f = Sigma_size;
    for (i=Sigma_size-1; i>=0; --i) {
      sum += Sigma[i] * Sigma[i];
      if (sum < epsqr)
	--cutoff.f;
      else
	break;
    }
	 
    // Extract the needed columns from U0 and write out to a file
    ptr_file = fopen(Umat,"wb");
    fwrite(&cutoff.f, sizeof(int), 1, ptr_file);
    fclose(ptr_file);
  }
  else {
    ptr_file = fopen(Umat, "rb");
    if (ptr_file == NULL)
      printf("Can't read cutoff.f!\n");
    READ_CHECK( fread(&cutoff.f, sizeof(int), 1, ptr_file), 1 );
    fclose(ptr_file);
  }
     
  //*Ucomp = ((double)cutoff.f)/Sigma_size;
  printf("U compression rate: %.4f\n", ((double)cutoff.f)/Sigma_size);

  double trancatedSize = dofn3_f * cutoff.f;
  ptr_file = fopen(Umat, "ab");

  // zero file
  //double *UU = calloc( trancatedSize, sizeof(double) );
  fwrite(U0, sizeof(double), trancatedSize, ptr_file);
  //fwrite(UU, sizeof(double), trancatedSize, ptr_file);
  fclose(ptr_file);
	
  /*
    char sinfile[50];
    sprintf(sinfile, "FieldSign%dA%.2f.out", n, AnisoParameters->elasA);
    ptr_file = fopen(sinfile, "w");
    fprintf(ptr_file, "%d\n", cutoff.f);

    for (i=0; i<Sigma_size; i++)
    fprintf(ptr_file, "%e\n", Sigma[i]);
    fclose(ptr_file);
  */
  free(Sigma); Sigma = NULL;
	
	
  /****
   * Compute the SVD of K_thin
   ****/
	
  // Form K_thin using all of 316 M2L operators stored in Kcell
  count = 0;
  for (j=0;j<dofn3_s;++j) {		
    for (i=0;i<316;i++) {
      for (l=0;l<dofn3_f;l++) {
	K0[count] = Kcell[i][l+j*dofn3_f];
	count++;
      }
    }
  }
	
  // save = "S"; nosave = "N"
  double *U1 = NULL;
  int rows_f = 316*dofn3_f;
  Sigma_size = dofn3_s;	
  Sigma      = (double *)malloc(Sigma_size * sizeof(double));
  VT0        = (double *)malloc(dofn3_s * dofn3_s * sizeof(double));

  // Compute the SVD of the K_thin
  dgesvd_(nosave,save,&rows_f,&dofn3_s,K0,&rows_f,Sigma,U1,&nosavedim,VT0,&dofn3_s,
	  work,&lwork,&info);

  if (callTime == 0) {
	 
    // Determine the number of singular values to keep. We use
    // epsilon for this.
	 
    sumsqr = 0;
    for (i=Sigma_size-1; i>=0; --i) 
      sumsqr += Sigma[i] * Sigma[i];

    sum = 0, epsqr = sumsqr * epsilon * epsilon;
    cutoff.s = Sigma_size;
    for (i=Sigma_size-1; i>=0; --i) {
      sum += Sigma[i] * Sigma[i];
      if (sum < epsqr)
	--cutoff.s;
      else
	break;
    }
	
    // Extract trancated VT[cutoff.s][dofn3_s] from
    // VT0[dofn3_s][dofn3_s] and write out to a file
    
    ptr_file = fopen(Vmat,"wb");
    fwrite(&cutoff.s, sizeof(int), 1, ptr_file);
    fclose(ptr_file);
  }
  else {
    ptr_file = fopen(Vmat, "rb");
    if (ptr_file == NULL)
      printf("Can't read cutoff.s!\n");
    READ_CHECK( fread(&cutoff.s, sizeof(int), 1, ptr_file), 1 );
    fclose(ptr_file);
  }
     
  //*Vcomp = ((double)cutoff.s)/Sigma_size;
  printf("V compression rate: %.4f\n", ((double)cutoff.s)/Sigma_size);

  ptr_file = fopen(Vmat, "ab");

  //double *VV = calloc( dofn3_s*cutoff.s, sizeof(double) );
  //fwrite(VV, sizeof(double), dofn3_s*cutoff.s, ptr_file);
      

  for (j=0;j<dofn3_s;j++) { // column
    count1 = j*dofn3_s;
    fwrite(VT0+count1, sizeof(double), cutoff.s, ptr_file);
  }

  fclose(ptr_file);

  /*
    sprintf(sinfile, "SourceSign%dA%.2f.out", n, AnisoParameters->elasA);
    ptr_file = fopen(sinfile, "w");
    fprintf(ptr_file, "%d\n", cutoff.s);
    for (i=0; i<Sigma_size; i++)
    fprintf(ptr_file, "%e\n", Sigma[i]);
    fclose(ptr_file);
  */
      
  /** Computing the compressed kernel using the orthonormal basis U and VT.
   **/
  char *transa, *transb;
  int cutoff2 = cutoff.f * cutoff.s;	
  double alpha=1, beta=0;
  double *Ecell, *KV;
  Ecell = (double *)malloc(cutoff2 * sizeof(double));	
  KV    = (double *)malloc(dofn3_f * cutoff.s * sizeof(double));
	
  ptr_file = fopen(Kmat,"ab");
  for (i=0;i<316;i++) {       

    /* Compute K V:
     * K  is dofn3_f  x dofn3_s
     * VT is cutoff.s x dofn3_s
     * V  is dofn3_s  x cutoff.s 
     * KV is dofn3_f  x cutoff.s
     * (Notice that VT is a submatrix of VT0)
     */
    transa = "n";
    transb = "t";
    dgemm_(transa, transb, &dofn3_f, &cutoff.s, &dofn3_s, &alpha, 
	   Kcell[i], &dofn3_f, VT0, &dofn3_s, &beta, KV, &dofn3_f);
		
    /* Compute U^T K V:
     * KV is dofn3_f x cutoff.s
     * U  is dofn3_f x cutoff.f
     * UT is cutoff.f x dofn3_f
     * U^T K V is cutoff.f x cutoff.s
     * (Notice that U is a submatrix of U0)
     */
    transa = "t";
    transb = "n";
    dgemm_(transa, transb, &cutoff.f, &cutoff.s, &dofn3_f, &alpha, 
	   U0, &dofn3_f, KV, &dofn3_f, &beta, Ecell, &cutoff.f);
		
    fwrite(Ecell, sizeof(double), cutoff2, ptr_file);
  }
  fclose(ptr_file);

  free(K0);
  free(VT0);
  free(U0);
  free(Sigma);
  free(ker);
  free(work);
  free(KV);
  free(Ecell);

  for (z = 0; z < 316; ++z)
    free(Kcell[z]);
     
}


/*
 * Function: ComputeKernelUnif
 * ------------------------------------------------------------------
 * Computes the kernel for 316(2n-1)^3 interactions between Uniform
 * Grid nodes.
 */
void ComputeKernelUnif(int n, int2 dof, kfun_t kfun, char *Kmat,
		       double alphaAdjust, double len) {
      
  int i, k1, k2, k3, l1, l2, l3;
  vec3 scenter;
    
  int dof2 = dof.s * dof.f;
  double nodes[n], kernel[dof2];
  vec3 fieldpos, sourcepos;
  GridPos1d(alphaAdjust, len, n, UNIF, nodes);

  int vecSize = 2*n-1, reducedMatSize = pow(vecSize, 3);
  int M2LSize = dof2 * reducedMatSize;
  double *MatM2L  = fftw_alloc_real(     M2LSize );
  double *freqMat = fftw_alloc_real( 316*M2LSize );
  fftw_plan p[316*dof2];
  
  // Compute the kernel values for interactions with all 316 cells
  int countM2L, countPlan, count, count1;
  int shift1, shift2, shiftGlo, shiftLoc;
  int reducedMatSizeDofs = reducedMatSize *dof.s;
  int f, s;

  countM2L  = 0;
  countPlan = 0;
  for (k1=-3;k1<4;k1++) {
    scenter.x = (double)k1*len;
    for (k2=-3;k2<4;k2++) {
      scenter.y = (double)k2*len;
      for (k3=-3;k3<4;k3++) {
	scenter.z = (double)k3*len;

	vec3 r = {k1, k2, k3};
	if ( Is_well_separated(r, 1.0) ) {

	  count=0;
	  for (l1=0; l1<vecSize; l1++) {
	    GetPosition(n, l1, &fieldpos.x, &sourcepos.x, nodes);
	    sourcepos.x += scenter.x;
	    for (l2=0; l2<vecSize; l2++) {
	      GetPosition(n, l2, &fieldpos.y, &sourcepos.y, nodes);
	      sourcepos.y += scenter.y;
	      for (l3=0; l3<vecSize; l3++, count++) {
		GetPosition(n, l3, &fieldpos.z, &sourcepos.z, nodes);
		sourcepos.z += scenter.z;
		
		kfun(&fieldpos, &sourcepos, kernel);

		// kernel[f x s] in column major
		// MatM2L[(2n-1)^3 x s x f] in column major
		count1 = 0;
		shift1 = count;
		for (s=0; s < dof.s; s++) {
		  shift2 = 0;
		  for(f=0; f < dof.f; f++) {
		    MatM2L[shift1 + shift2] =  kernel[count1++];
		    shift2 += reducedMatSizeDofs;
		  }
		  shift1 += reducedMatSize;
		}
	      }
	    }
	  }	      
	      
	  // FFT
	  shiftGlo = countM2L *M2LSize;

	  for (i=0, shiftLoc=0; i<dof2; i++) {
	    //rfftw_one(fft_r2c, MatM2L + shiftLoc, freqMat + shiftGlo + shiftLoc);

	    p[ countPlan ] = fftw_plan_r2r_1d(reducedMatSize, MatM2L + shiftLoc, freqMat + shiftGlo + shiftLoc, FFTW_R2HC, FFTW_FLAG);
	    fftw_execute( p[ countPlan ] );
	    countPlan++;
		
	    shiftLoc += reducedMatSize;
	  }

	  countM2L++;
	  //printf("countM2L: %d\n", countM2L++);
	}
      }
    }
  }


  FILE *ptr_file;
  ptr_file = fopen(Kmat, "ab");
  fwrite(freqMat, sizeof(double), 316 *M2LSize, ptr_file);
  fclose(ptr_file);

      
  fftw_free(MatM2L);
  fftw_free(freqMat);
  assert(countPlan == 316*dof2);
  for (i=0; i<countPlan; i++)
    fftw_destroy_plan( p[i] );
}

 
/*
 * Function: NewNode
 * ------------------------------------------------------------------
 * Dynamically allocates space for a new node A in the octree and
 * initializes the quantities in the node.
 *
 */
nodeT* NewNode(vec3 center, double L, int n) {
	
  // Initialization
  nodeT *A = (nodeT *)malloc(sizeof(nodeT));
	
  // Initializes the child, neighbor, interaction, and parent nodes to NULL
  int i;	
  for (i=0;i<8;i++)
    A->leaves[i] = NULL;
  A->parent = NULL;
	
  A->fieldval   = NULL;
  A->sourceval  = NULL;
  A->proxysval  = NULL;
  A->sourcefre  = NULL;

  A->fieldlist  = NULL;
  A->sourcelist = NULL;
  A->center     = center;
  A->length     = L;
  A->Nf         = 0;
  A->Ns         = 0;
  A->ineigh     = 0;
  A->iinter     = 0;

  return A;
}

/*
 * Function: BuildFMMHierarchy
 * ------------------------------------------------------------------
 * Builds the FMM hierarchy with l levels.
 *
 */
void BuildFMMHierarchy(nodeT *A, int level, int n) {
   
  double L, halfL, quarterL;
  vec3 center, left, right;
	
  if (level > 0) {
    // Compute the length and center coordinate of the children cells
    L        = A->length;
    halfL    = 0.5*L;
    quarterL = 0.25*L;
    center   = A->center;
    left.x   = center.x - quarterL;
    left.y   = center.y - quarterL;
    left.z   = center.z - quarterL;
    right.x  = center.x + quarterL;
    right.y  = center.y + quarterL;
    right.z  = center.z + quarterL;    
		
    // Add children nodes to node A
    int i;
    for (i=0;i<8;i++) {
      // Determine the center of the child cell
      if (i<4) {
	center.x = left.x;
				
	if (i<2)
	  center.y = left.y;
	else
	  center.y = right.y;
				
      } else {
	center.x = right.x;
				
	if (i<6)
	  center.y = left.y;
	else
	  center.y = right.y;
      }
			
      if (i%2 == 0)
	center.z = left.z;
      else
	center.z = right.z;
			
      // Create the child cell
      A->leaves[i] = NewNode(center,halfL,n);
      A->leaves[i]->parent = A;
    }
		
    // Recursively build octree if there is a subsequent level
    for (i=0;i<8;i++)
      BuildFMMHierarchy(A->leaves[i],level-1,n);
  }
}

	
/*
 * Function: DistributeFieldPoints
 * -------------------------------------------------------------------
 * Distributes all of the field points to the appropriate cell at the 
 * bottom of the octree.
 */
void DistributeFieldPoints(nodeT *A, vec3 *field, int *fieldlist, 
			   int levels) {
  int i, j, k, m, z;
  int Nf = A->Nf;
  vec3 point, center;
  int *fieldcell[8], *F;
  for (z = 0; z < 8; z++)
    fieldcell[z] = (int *)malloc(Nf * sizeof(int));
	
  if (levels > 0) {
    // Obtain the center of the cell
    center = A->center;
		
    // Distribute the field points    
    if (Nf != 0) {
      // Determine which child cell each field and source point belongs to
      for (i=0;i<Nf;i++) {
	k = fieldlist[i];
	point = field[k];
				
	// Determine which cell the point is in
	if (point.x < center.x) {
	  if (point.y < center.y) {
	    if (point.z < center.z)
	      j = 0;
	    else
	      j = 1;
	  } else {
	    if (point.z < center.z)
	      j = 2;
	    else
	      j = 3;
	  }
	} else {
	  if (point.y < center.y) {
	    if (point.z < center.z)
	      j = 4;
	    else
	      j = 5;
	  } else {
	    if (point.z < center.z)
	      j = 6;
	    else
	      j = 7;
	  }
	}
				
	// Add the field point to the list for child cell j
	m = A->leaves[j]->Nf;
	fieldcell[j][m] = k;
	A->leaves[j]->Nf++;
      }
			
      // Recursively distribute the points
      for (j=0;j<8;j++) {
	F = fieldcell[j];
	DistributeFieldPoints(A->leaves[j],field,F,levels-1);
      }
    }
  } else if (levels == 0) {
    Nf = A->Nf;
    A->fieldlist = (int *)malloc(Nf*sizeof(int));
    F = A->fieldlist;
    for (i=0;i<Nf;i++) 
      F[i] = fieldlist[i];
  }   
  for (z = 0; z < 8; z++)
    free(fieldcell[z]); 
}

/*
 * Function: DistributeSources
 * -------------------------------------------------------------------
 * Distributes all of the sources to the appropriate cell at the 
 * bottom of the octree.
 */
void DistributeSources(nodeT *A, vec3 *source, int *sourcelist, 
		       int levels) {
  int i, j, k, m, z;
  int Ns = A->Ns;
  vec3 point, center;
  int *sourcecell[8], *S;
  for (z = 0; z < 8; z++)
    sourcecell[z] = (int *)malloc(Ns * sizeof(int));
	
  if (levels > 0) {
    // Obtain the center of the cell
    center = A->center;
		
    // Distribute the field points
    if (Ns != 0) {      
      // Determine which child cell each field and source point belongs to
      for (i=0;i<Ns;i++) {
	k = sourcelist[i];
	point = source[k];
				
	// Determine which cell the point is in
	if (point.x < center.x) {
	  if (point.y < center.y) {
	    if (point.z < center.z)
	      j = 0;
	    else
	      j = 1;
	  } else {
	    if (point.z < center.z)
	      j = 2;
	    else
	      j = 3;
	  }
	} else {
	  if (point.y < center.y) {
	    if (point.z < center.z)
	      j = 4;
	    else
	      j = 5;
	  } else {
	    if (point.z < center.z)
	      j = 6;
	    else
	      j = 7;
	  }
	}

	//printf("j: %d\n", j);
	    
	// Add the field point to the list for child cell j
	m = A->leaves[j]->Ns;
	sourcecell[j][m] = k;
	A->leaves[j]->Ns++;
      }
			
      // Recursively distribute the points
      for (j=0;j<8;j++) {
	S = sourcecell[j];
	DistributeSources(A->leaves[j],source,S,levels-1);
      }
    }
  } else if (levels == 0) {
    Ns = A->Ns;
    A->sourcelist = (int *)malloc(Ns*sizeof(int));
    S = A->sourcelist;
    for (i=0;i<Ns;i++) 
      S[i] = sourcelist[i];
  }
  for (z = 0; z < 8; z++)
    free(sourcecell[z]);
}


/*
 * Function: InteractionList
 * -------------------------------------------------------------------
 * For each node at each level of the octree, the interaction list
 * is determined.
 */
void InteractionList(nodeT *A, int levels) {
  int i, j, k, ineigh, iinter, nneigh;
  nodeT *B, *C;
  vec3 center1, center2, cshift, diff;
  double L;
	
  assert(A->Nf > 0);
	
  if (levels > 0) {
    nneigh = A->ineigh;
		
    /* Sets the cutoff between near and far to be L
     * (this is equivalent to a one cell buffer of the child)
     */
    L = A->length / 2.0;
		
    /* 
     * Finds all neighbors that are too close for the far field 
     * approximation and stores them in the neighbors array - 
     * Also finds the neighboring cells that are sufficiently 
     * far away and stores them in interaction array
     */
    // loop over all the neighbors of A
    for (i=0;i<nneigh;i++) {
      B = A->neighbors[i]; 
      cshift = A->cshiftneigh[i];

      // loop over all the children of B
      for (j=0;j<8;j++) {
	if (B->leaves[j]->Ns != 0) { 
					
	  center1 = B->leaves[j]->center;
	  center1.x += cshift.x;
	  center1.y += cshift.y;
	  center1.z += cshift.z;

	  // loop over all the children of A
	  for (k=0;k<8;k++) {
	    C = A->leaves[k];
						
	    if (C->Nf != 0) {
							
	      ineigh  = C->ineigh;
	      iinter  = C->iinter;
	      center2 = C->center;
							
	      diff.x  = center1.x - center2.x;
	      diff.y  = center1.y - center2.y;
	      diff.z  = center1.z - center2.z;
							
	      if ( Is_well_separated(diff, L) ) {
		C->interaction[iinter] = B->leaves[j];
		C->cshiftinter[iinter] = cshift;
		C->iinter++;

	      } else { // This is a neighbor
		C->neighbors[ineigh] = B->leaves[j]; 
		C->cshiftneigh[ineigh] = cshift;
		C->ineigh++;
	      }
	    }
						
	  }
	}
      }
    }
		
    // recursively build the interaction lists
    if (A->leaves[0] != NULL) {
      for (k=0;k<8;k++) {
	if (A->leaves[k]->Nf != 0) {
	  InteractionList(A->leaves[k],levels-1);
	}
      }
    }
  }
}
 
 
/*
 * Function: Upward Pass
 * Recursive call M2M - SVD compression Step 8718 - 3a
 * Loop over the level - Keep that - Tkz : info about Chebychev's coeffs.
 * Independent of the Kernel
 *
 * Input:
 *
 *   	   A: Pointer to the pointer of a box in the hierarchy tree so *A is 
 *            the pointer to the box (see below for variablesused in the structure,
 *            refer to file bbfmm.h for the nodeT structure)
 *
 * (*A)->leaves[j]->sourceval: Source values of children boxes (refer to W on the
 *                             right hand side in step 2 at the bottom of page 8716)
 *
 * (*A)->sourcelist: Pointer to the array of global indices of sources in the box, used
 *                   with 'segment' below to find the position of segments in the box
 *
 * (*A)->Ns: Number of sources in the box
 *
 * (*A)->center: Position of the box center 
 *
 * (*A)->length: Length of the box
 *
 * (*A)->leaves: Pointer to the children boxes, used to determine whether A is a leaf by
 *               checking (*A)->leaves[0] == NULL, and recursivly go up the tree
 *
 * (*A)->leaves[i]->Ns: Number of source points in the child box, used to determine whether
 *                      the it is an empty box to avoid computation by checking
 *                      (*A)->leaves[i]->Ns != 0
 *
 *   segment: Pointer to the segment array where every segment contains
 *            two end points (to be used in P2M)
 *
 *  midpoint: Pointer to the array of midpoints of all the segments 
 *            (to be used in P2M)
 *
 *  Cweights: Pointer to the array of Chebyshev interpolation coefficients, 
 *  	      transfering the weights at Chebyshev nodes between children and parent
 *            (Refer to Sn on the right hand side in step 2 and step 4 
 *            at bottom of page 8715)
 *
 *       Tkz: Pointer to the n x n array (n is the number of chebyshev nodes) 
 *            containing the evaluation of n Chebyshev polynomials 
 *            (T_0 to T_{n-1}) at n chebyshev nodes (of T_n)
 *
 *      burg: Pointer to the array of Burger's vectors (to be used in P2M)
 *
 *        VT: Pointer to the transpose of the truncated right singular
 *            vectors (refer to S_r^k in step 0 at page 8718)
 *
 *    cutoff: cutoff.f is the cutoff of the left singular vectors U and 
 *            cutoff.s is the cutoff of the right singular vectors V
 *    
 *	   n: number of Chebyshev nodes of every dimension in every box 
 *
 *       dof: dof.f is the degree of the field and dof.s is the degree
 *            of the source
 *
 *    numgau: Number of Gaussian points used in line integral (to be used 
 *            in P2M)
 *
 *    gpoint: Gaussian points in [-1, 1] (to be used in P2M)
 *
 *   gweight: Gaussian quadrature weights (to be used in P2M)
 *
 *
 * Output:
 * (*A)->sourceval: Source value at the Chebyshev nodes (refer to W on the left hand side
 *                  in step 2 at the bottom of page 8716 and w^l on the right hand side
 *                  in step 3a at page 8718)
 *
 * ---------------------------------------------------------------------
 * Gathers the sources from the children cells and determines the strength
 * of pseudo-charges located at the Chebyshev nodes of the parent cell.
 * (upward pass of BBFMM)
 */
void UpwardPass(nodeT **A, FMMSrc fmm_src, double *Cweights,
		double *Tkz, int n, int dof, double alpha, int grid_type) {

  
  int i, l;
  int dofn3 = dof * n*n*n;

	
  /* First gather the sources for all children cells and
   * then gather for the parent cell - otherwise map all
   * sources to Chebyshev nodes */
  if ((*A)->leaves[0] != NULL) {  // Going up the tree M2M

    // Allocate memory for the source values and initialize
    (*A)->sourceval = (double *)calloc(dofn3, sizeof(double));
	  
    // Determine which children cells contain sources
    for (i=0;i<8;i++) {
      if ((*A)->leaves[i]->Ns != 0) { 


	/*
	  #ifdef LINEINT        
	  // Recursively gather sources
	  UpwardPass(&((*A)->leaves[i]),segment,midpoint,Cweights,Tkz,burg,
	  n,dof,gpoint,gweight, alpha, grid_type);
	                
	  #elif TENSOR
	*/
	UpwardPass(&((*A)->leaves[i]), fmm_src, Cweights, Tkz,
		   n, dof, alpha, grid_type);
	

	double *Schild;      
	Schild = (*A)->leaves[i]->sourceval;
			
	// Manually set the vector from child to parent
	int  idx[3] = {i/4%2, i/2%2, i%2};
	double r[3];
	for (l=0; l<3; l++)
	  if (idx[l] == 0)
	    r[l] = -1;
	  else
	    r[l] =  1;

	double SS[dofn3];
	Moment2Moment(n, r, Schild, SS, dof, Cweights);


	// the child's source value will not be used any more
	if (grid_type == UNIF) {
	  free(Schild), Schild = NULL;
	  (*A)->leaves[i]->sourceval = NULL;
	}

	(*A)->sourceval = Add((*A)->sourceval, SS, dofn3);
      }
    }
	  
  } else { // leaf boxes

    if (fmm_src.src_t == SEG) // optimized P2M for segments
      P2M_SEG(A, fmm_src, n, dof, Tkz, alpha, grid_type);
    else 
      P2M_PNT(A, fmm_src, n, dof, Tkz, grid_type);	
  }
	
  if (grid_type == UNIF) {
	
    // pad the source value
    int padSize = (int)round(pow(2*n-1, 3));
    int halfSize = padSize/2;
    int l3Size = n, l1Size = (int)round(pow(l3Size, 2));
    int l2Pad = 2*n-1, l1Pad = (int)round(pow(l2Pad, 2));
    int shift, count;
    int s, l1, l2, l3, n3 = n*n*n;

    double *x = (*A)->sourceval;
    //double *padx = calloc(padSize * dof, sizeof(double));
    //(*A)->sourcefre = malloc(padSize * dof *sizeof(double));
    double *padx    = fftw_alloc_real( padSize * dof );
    (*A)->sourcefre = fftw_alloc_real( padSize * dof );

    Zero(padx, padSize * dof);
    for (i=0, count=0; i<n3; i++) {
      l3 = i % l3Size;
      l1 = i / l1Size;
      l2 = i / l3Size % l3Size;
      shift = halfSize+l1*l1Pad+l2*l2Pad+l3;
	    
      for (s=0; s<dof; s++) {
	padx[shift] = x[count++];
	shift += padSize; 
      }
    }

    fftw_plan p[dof];
    for (s=0, shift=0; s<dof; s++) {
      //rfftw_one(fft_r2c, padx + shift, (*A)->sourcefre + shift);
      p[s] = fftw_plan_r2r_1d(padSize, padx + shift, (*A)->sourcefre + shift, FFTW_R2HC, FFTW_FLAG);
      fftw_execute( p[s] );
      shift += padSize;
    }
	
    x=NULL;
    fftw_free(padx), padx=NULL;
    for (s=0; s<dof; s++)
      fftw_destroy_plan(p[s]);
  }

}

	    
/* Function: M2M 
 * Input:
 *       n: Number of Chebyshev nodes
 *       r: (r[0], r[1], r[2]) is the vector from the child
 *          box to the parent box
 *  Schild: Source values of the child
 *
 * Output:
 *  SS: Source values of the parent
 *
 */

void Moment2Moment(int n, double *r, double *Schild, double *SS, 
		   int dof, double *Cweights) {
     

  int l, l1, l2, l3, l4;
  int count1, count2, count3, wstart;
	
  int incr  = 1;
  int n2    = n*n;                       // n2 = n^2
  int n3    = n*n*n;                     // n3 = n^3	
  int dofn  = dof * n;
  int dofn2 = dof * n2;
  int dofn3 = dof * n3;
	
  double Sy[dofn3], Sz[dofn3]; 
    

  // Recover the index of the child box
  int idx[3] = {0, 0, 0};
  if (  r[0] > 0 )
    idx[0] = 1;
  if (  r[1] > 0 )
    idx[1] = 1;
  if (  r[2] > 0 )
    idx[2] = 1; 

  // Gather the children source along the z-component
  if (idx[2] == 0)
    wstart =  0;
  else
    wstart =  n2;
			
  l = 0;
  for (l1=0;l1<n2;l1++) {
    count1 = l1*n;
    for (l3=0;l3<n;l3++) {
      count2 = wstart + l3*n;
      for (l4=0;l4<dof;l4++) {
	count3 = dof*count1 + l4;
	Sz[l]  = ddot_(&n, &Schild[count3], &dof, &Cweights[count2], &incr);
	l++;
      }
    }
  }
		
  // Gather the children sources along the y-component
  if (idx[1] == 0)
    wstart =  0;
  else
    wstart =  n2;
			
  l = 0;
  for (l1=0;l1<n;l1++) {
    for (l2=0;l2<n;l2++) {
      count2 = wstart + l2*n;
      for (l3=0;l3<n;l3++) {
	count1 = l1*n2 + l3;
	for (l4=0;l4<dof;l4++) {
	  count3 = dof*count1 + l4;
	  Sy[l]  = ddot_(&n,&Sz[count3],&dofn,&Cweights[count2],&incr);
	  l++;
	}
      }
    }
  }
		
  /* Gather the children sources along the x-component and determine
     the parent sources */
  l = 0;
  if (idx[0] == 0)
    wstart =  0;
  else
    wstart =  n2;
			
  for (l1=0;l1<n;l1++) {
    count2 = wstart + l1*n;
    for (l2=0;l2<n;l2++) {
      for (l3=0;l3<n;l3++) {
	count1 = l2*n + l3;
	for (l4=0;l4<dof;l4++) {
	  count3  = dof*count1 + l4;
	  SS[l] = ddot_(&n, &Sy[count3], &dofn2, &Cweights[count2], &incr);
	  l++;
	}
      }
    }
  }
}


void P2M_PNT(nodeT **A, FMMSrc fmm_src, int n, int dof, double *Tkz, gridT grid_type) {
  
  int Ns = (*A)->Ns; // Number of source points in a leaf box
  int j, k, l1, l2, l3, l4, *sourcelist = (*A)->sourcelist;
  double ihalfL = 2./(*A)->length;
  vec3 scenter = (*A)->center;
  vec3 *sourcet = malloc(  Ns * sizeof(vec3));
  vec3 *Ss      = malloc(n*Ns * sizeof(vec3));
  
  vec3 *source = fmm_src.source;
  double *q = fmm_src.q;
    
  // Map all of the source points to the box ([-1 1])^3
  for (j=0;j<Ns;j++) {
    k = sourcelist[j];
    sourcet[j].x = ihalfL*(source[k].x - scenter.x);
    sourcet[j].y = ihalfL*(source[k].y - scenter.y);
    sourcet[j].z = ihalfL*(source[k].z - scenter.z);
  }
	
  // Compute Ss, the mapping function for the sources
  ComputeSn(sourcet, Ns, Ss, n, Tkz, grid_type);
        
  // Compute the source values
  int dofn3 = dof * n*n*n;
  (*A)->sourceval = (double *)calloc(dofn3,sizeof(double));
  double *S = (*A)->sourceval;

  int x_idx, y_idx, z_idx, lhs_idx, rhs_idx;
  double weight;

  lhs_idx = 0;
  for (l1=0;l1<n;l1++) {
    x_idx = l1*Ns;
    for (l2=0;l2<n;l2++) {
      y_idx = l2*Ns;
      for (l3=0;l3<n;l3++) {
	z_idx = l3*Ns;

	for (j=0; j<Ns; j++) {
	  weight = Ss[x_idx+j].x *Ss[y_idx+j].y *Ss[z_idx+j].z;
	  rhs_idx = sourcelist[j] * dof;
		  
	  for (l4=0; l4<dof; l4++)
	    S[lhs_idx+l4] += weight * q[rhs_idx++];

	}
	lhs_idx += dof;
		
      }
    }
  }
}


/* Function: P2M
 *   Input:
 *
 *   	   A: Pointer to the pointer of a box in the hierarchy tree so *A is 
 *            the pointer to the box (see below for variables used in the structure,
 *            refer to file bbfmm.h for the nodeT structure)
 *
 * (*A)->sourcelist: Pointer to the array of global indices of sources in the box,
 *                   used with 'segment' below to find the position of segments in the box
 *
 * (*A)->Ns: Number of sources in the box
 *
 * (*A)->center: Position of the box center 
 *
 * (*A)->length: Length of the box
 *
 *   segment: Pointer to the segment array where every segment contains
 *            two end points
 *
 *      burg: Pointer to the array of Burger's vectors (to be used in P2M)
 *
 *  midpoint: Pointer to the array of midpoints of all the segments 
 *
 *         n: Number of Chebyshev nodes of every dimension in every box 
 *
 *       Tkz: Pointer to the n x n array (n is the number of chebyshev nodes) 
 *            containing the evaluation of n Chebyshev polynomials 
 *            (T_0 to T_{n-1}) at n chebyshev nodes (of T_n)
 *
 *        VT: Pointer to the transpose of the truncated right singular
 *            vectors (refer to S_r^k in step 0 at page 8718)
 *
 *       dof: dof.f is the degree of the field and dof.s is the degree
 *            of the source
 *
 *    numgau: Number of Gaussian points used in line integral (to be used 
 *            in P2M)
 *
 *    gpoint: Gaussian points in [-1, 1] (to be used in P2M)
 *
 *   gweight: Gaussian quadrature weights (to be used in P2M)
 *
 * Output:
 *
 * (*A)->sourceval: Source value at the Chebyshev nodes (refer to W on the
 * left hand side
 *                  in step 1 at the bottom of page 8716 and w^l on the
 * right hand side in
 *                  step 3a at page 8718)
 *	
 */

void P2M_SEG(nodeT **A, FMMSrc fmm_src, int n, int dof, double *Tkz,
	     double alpha, int grid_type) {

  int     nGauss = fmm_src.nGauss;
  double *burg   = fmm_src.burger;
  
  int  Ns         = (*A)->Ns;             // Number of source points
  int  dofn3      = dof * n*n*n;
  int  numSource  = Ns  * nGauss;
  int *sourcelist = (*A)->sourcelist; 	  // Indices of sources in cell
  vec3 scenter    = (*A)->center;    	  // Center of source cell
  vec3 *Ss        = (vec3 *)malloc(n * numSource * sizeof(vec3));	
	
    
  // Adjust box size to accomodate outside segments
  double L = AdjustBoxSize((*A)->length, alpha);
    
  // Compute source points from every segment in the unit box
  double *xi = (double *) malloc(3 * Ns * sizeof(double)); 
  vec3 *sourcet = (vec3 *)malloc(numSource * sizeof(vec3));
  LineIntegral(fmm_src, sourcelist, Ns, L, scenter, sourcet, xi);
    
  // Compute Ss, the mapping function for the sources
  ComputeSn(sourcet, numSource, Ss, n, Tkz, grid_type);
    
  // Compute the source values
  (*A)->sourceval = (double *)calloc(dofn3,sizeof(double));
  double *S = (*A)->sourceval;
    
  /* How is dof used?
     q stores many right hand sides (columns). That is we are computing
     a matrix-vector product A q, where A contains the kernel and q is a matrix
     with dof columns. 
     The storage in q is row-major, that is the first row is q[0], ..., q[dof-1].
     The second row is q[dof], ..., q[2*dof-1]. And so on.
  */

  int j, k, l1, l2, l3, l5, indx;
  double sum;

  int lhs_idx = 0, idof;
  for (l1=0;l1<n;l1++) {
    for (l2=0;l2<n;l2++) {
      for (l3=0;l3<n;l3++) {
	      
	for (j=0;j<Ns;j++) {
	  sum = 0;
	  for (l5=0;l5<nGauss;l5++) {
	    /* Row major storage of data in q
	     * row = sourcelist[j]
	     * column = l4
	     */
	    // page 8716 Eq. 1
	    indx = j*nGauss + l5;
	    sum += GAUSSW[l5] * Ss[l1*numSource+indx].x *
	      Ss[l2*numSource+indx].y * Ss[l3*numSource+indx].z;
			    
	  }
	  // Note that 'xi' is local of the box while 'burg' is global
	  k = sourcelist[j];
		
	  /* Row major storage in S as well.
	   * S[0], ..., S[dof-1] corresponds to multipole
	   * coefficient 0 for all dofs.
	   */
	  for (idof=0; idof<dof; idof++)
	    S[lhs_idx+idof] += sum * burg[k*3 + idof%3] * xi[j*3 + idof/3];

	}
	lhs_idx += dof;
      }	
    }
  }
  free(xi);
  free(Ss);
  free(sourcet);

  /*
  printf("\n\n");
  printf("Ns: %d\n", Ns);
  print_array((*A)->sourceval, dofn3, "Source weights");
  */
}

      
/*
 * return: L and phi
 */
void AddPBC( int lpbc, double *M, double **L, double *phi, int Nf,
	     int n, int2 dof, double boxLen, double alpha, double *Tkz,
	     kernel_t kernel, int grid_type ) {

  int n3  = n*n*n;
  int n3f = n3*dof.f;
  (*L) = calloc( n3f, sizeof(double) ); // initialize PBC stress
      
  if (lpbc > 1) { // Add PBC stress

    printf("Start PBC ...\n");
     
    int n3s = n3*dof.s;
    double *Unif, *Cheb; // values on the uniform and Chebyshev grids
    if (grid_type == UNIF) {
      Unif = M;
      Cheb = calloc( n3s, sizeof(double) );
      Anterpolate(Unif, UNIF, Cheb, CHEB, n, dof.s, Tkz);
      M    = Cheb; // source values on the Chebyshev grids
    }
    
    // compute PBC potential on the Chebyshev grids
    DoCorrectionTable(M, (*L), n, dof, boxLen, alpha, lpbc, kernel, Tkz);
    //print_array((*L), n3f, "FMM PBC on Chebyshev");
    
    // compute PBC mean from the values on the Chebyshev grids
    double *PBCmean = MeanStressFmm((*L), n, dof.f, Tkz);    
    int Nff = Nf*dof.f;
    phi = Subtract( phi, Nff, PBCmean, dof.f ); // subtract mean
    free(PBCmean);

    if (grid_type == UNIF) {
      free(Cheb);  // source values on the Chebyshev grids
      Cheb = (*L); // field  values on the Chebyshev grids
      Unif = calloc( n3f, sizeof(double) );
      Interpolate(Cheb, CHEB, Unif, UNIF, n, dof.f, Tkz);
      free(Cheb);
      (*L) = Unif;
    }
	    
  }

}

      
/*
 * Function: FMMInteraction
 *   Input:
 *
 *   	   A: Pointer to the pointer of a box in the hierarchy tree so *A is 
 *            the pointer to the box (see below for variables used in the structure,
 *            refer to file bbfmm.h for the nodeT structure)
 *
 * (*A)->interaction[i]->sourceval: Pointer to the array of compressed values on source
 *                                  Chebyshev nodes (refer to w on the right hand side in
 *                                  step 3b at page 8718)
 *
 * (*A)->interaction[i]->center: Position of the interaction box center
 *
 * (*A)->leaves: Pointer to leaf boxes, used to recursively go over the tree
 *
 * (*A)->leaves[i]->Nf: Number of field points in the leaf box, used to check whether the box is empty
 *
 * (*A)->cshiftinter: used in PBC with 'Ktable' below for the indice of the M2L operator
 *
 * (*A)->center: Position of the box center 
 *
 * (*A)->length: Length of the box
 *
 * (*A)->iinter: Number of (nonempty) interaction boxes among the 189 candidates
 *
 * (*A)->interaction: Pointer to the array of pointers to the interaction boxes, so
 *                    (*A)->interaction[i] gives the pointer to the interaction box
 *
 *	   E: Pointer to the array of all the compressed M2L operators (refer to step 0
 *            at the page of 8718)
 *
 *    Ktable: Pointer to the array of indices of M2L operators, used to find a specific
 *            compressed M2L operator in 'E' above 
 *
 *         U: Pointer to the array of truncated left singular vectors (refer to the SVD
 *            of K_fat at page 8718)
 *
 *  Kweights: (refer to omega matrix in eqn. (8) at page 8717)
 *  
 *       dof: dof.f is the degree of the field and dof.s is the degree
 *            of the source
 *
 *    cutoff: cutoff.f is the cutoff of the left singular vectors U and 
 *            cutoff.s is the cutoff of the right singular vectors V
 *    
 *   Homogen: Order of the homogeneous kernel
 *
 *  Output:
 *  	(*A)->fieldval: Field values at the Chebyshev nodes (refer to the left hand side
 *                      in step 3c at page of 8719)
 *
 * ---------------------------------------------------------------------
 * At each level of the octree the interaction between well-separated cells of 
 * field and source Chebyshev nodes is computed using a SVD lookup table.
 */
void FMMInteraction(nodeT **A, double *E, int *Ktable, double *U, 
		    double *VT, double *Kweights, int n, int2 dof,
		    int2 cutoff, double homogen, int curTreeLevel, int
		    grid_type) {
  
  int lvl_shift;
  int Ksize;
  if ( !IsHomoKernel(homogen) ) {
    lvl_shift = (curTreeLevel>=2) ? curTreeLevel-2: 0;
    Ksize = 316*(2*n-1)*(2*n-1)*(2*n-1)*dof.s*dof.f;
     
  } else {
    lvl_shift = 0;
    Ksize = 316 * cutoff.f * cutoff.s;
  }
     
  int i, j, l;
  int n3       = n*n*n; 
  int dofn3_f  = dof.f * n3;
  int dofn3_s  = dof.s * n3;
  int cutoff_f = cutoff.f;
  int Usize    = dofn3_f * cutoff.f;
  int Vsize    = dofn3_s * cutoff.s;
  double L     = (*A)->length;    // Length of cell
  double iL    = 1.0/L;           // Inverse-length  
  double scale = pow(iL,homogen); // Scaling factor for M2L
  
  assert((*A)->Nf > 0); /* Cell cannot be empty. */
     
  double *productfre = NULL; // for uniform grids
  int matSizeDof;
  if (grid_type == CHEB)
    matSizeDof = cutoff_f;
  else if (grid_type == UNIF) {
    matSizeDof = (int)round(pow(2*n-1,3)) *dof.f;
    productfre = fftw_alloc_real( matSizeDof );
    Zero(productfre, matSizeDof);
  } else
    assert( false );
	
	
  double *FFCoeff = malloc(matSizeDof *sizeof(double));
	
  // Allocate memory for field values
  if ((*A)->parent == NULL) {
    assert( (*A)->fieldval != NULL );
    //print_array( (*A)->fieldval, n*n*n*dof->f, "root field value");
  } else {
    assert( (*A)->fieldval == NULL );
    (*A)->fieldval = (double *)calloc(dofn3_f, sizeof(double)); // initialize to zero
    assert( (*A)->fieldval != NULL );
  }
     
  vec3 fcenter = (*A)->center;   // Obtain the center of the field cell
  int ninter   = (*A)->iinter;   // Obtain the number of interaction cells
  double *F    = (*A)->fieldval; // Initialize pointer to the field values of A
  double *Pf   = calloc(cutoff_f, sizeof(double)); // Initialize to zero
	
  // Compute the field values due to all members of the interaction list
  for (i=0;i<ninter;i++) {	
		
    nodeT *B = (*A)->interaction[i]; 

    // Obtain the center of the source cell
    vec3 scenter = B->center;
    vec3 cshift  = (*A)->cshiftinter[i];
    scenter.x += cshift.x;
    scenter.y += cshift.y;
    scenter.z += cshift.z;
		
    // Note this is the normalized vector
    double R[3] = { iL*(scenter.x-fcenter.x),
		    iL*(scenter.y-fcenter.y),
		    iL*(scenter.z-fcenter.z) };

    // initialize to zeros
    for (j=0; j<matSizeDof; j++) {
      FFCoeff[j] = 0;
    }

    if (grid_type == CHEB) {
       
      Moment2Local(n, R, B->sourceval, FFCoeff, E + Ksize*lvl_shift, Ktable,
		   dof, cutoff, VT + Vsize*lvl_shift, Kweights, grid_type);	      
      Pf = Add(Pf, FFCoeff, cutoff_f);
    }
		
    else if (grid_type == UNIF) {
      Moment2Local(n, R, B->sourcefre, FFCoeff, E + Ksize*lvl_shift, Ktable, dof, cutoff, VT, Kweights, grid_type);
      productfre = Add(productfre, FFCoeff, matSizeDof);
    }
  } // end for ninter

  free(FFCoeff), FFCoeff=NULL;
     
  if (grid_type == CHEB) {

    int incr     =  1;
    char trans   = 'n';
    double alpha =  0;
    double F_m2l[n*n*n*dof.f];
    dgemv_(&trans, &dofn3_f, &cutoff_f, &scale, U + Usize*lvl_shift,
	   &dofn3_f, Pf, &incr, &alpha, F_m2l, &incr);

    
    // Adjust the field values by the appropriate weight
    l = 0;
    for (i=0;i<n3;i++) {
      double tmp = Kweights[i];
      for (j=0;j<dof.f;j++) {
	F_m2l[l] *= tmp;
	l++;
      }
    }

    F = Add(F, F_m2l, n*n*n*dof.f); // F of the root box has PBC stress
    
  } else if (grid_type == UNIF) { // uniform grids
	    
    int padSize = round(pow(2*n-1, 3));
    int l3Size = n, l1Size = (int)round(pow(l3Size, 2));
    int l2Pad = 2*n-1, l1Pad = (int)round(pow(l2Pad, 2));
    int shift, count;
    int f, l1, l2, l3;
	 
    double *res = fftw_alloc_real( padSize*dof.f );
    fftw_plan p[dof.f];
    for (f=0; f<dof.f; f++) {
      //rfftw_one(fft_c2r, productfre + f*padSize, res + f*padSize);
      p[f] = fftw_plan_r2r_1d(padSize, productfre + f*padSize, res + f*padSize, FFTW_HC2R, FFTW_FLAG);
      fftw_execute(p[f]);
    }

    for (f=0; f<dof.f; f++)
      fftw_destroy_plan(p[f]);

     
    fftw_free(productfre), productfre=NULL;
	    
    for (i=count=0; i<n3; i++) {
      l3 = i % l3Size;
      l1 = i / l1Size;
      l2 = i / l3Size % l3Size;
      shift = l1*l1Pad+l2*l2Pad+l3;
	     
      for (f=0; f<dof.f; f++, shift+=padSize)
	Pf[count++] = res[shift]/padSize;
	     
    }

    fftw_free(res), res = NULL;
     
    // F += FFT result
    for (j=0; j<cutoff_f; j++)
      F[j] += scale *Pf[j]; 
  }
     
  free(Pf);

  
  // Recursively compute the kernel interactions for all children cells
  // - 3c p. 8718
  // Go to the next level
  if ((*A)->leaves[0] != NULL) {
    for (i=0;i<8;i++) {
      if ((*A)->leaves[i]->Nf != 0) 
	FMMInteraction(&((*A)->leaves[i]), E, Ktable, U, VT,
		       Kweights, n, dof, cutoff, homogen,
		       curTreeLevel+1, grid_type);
    }
  }
}

/* Function: M2L
 * Input:
 * cell_mpCoeff:
 *
 *	   E: (Precomputed) Pointer to the array of all the compressed M2L operators 
 *	      (refer to step 0 at the page of 8718)
 *
 *    Ktable: (Precomputed) Pointer to the array of indices of M2L operators, used to 
 *            find a specific compressed M2L operator in 'E' above
 *
 * Output:
 *   FFCoeff:
 */
void Moment2Local(int n, double *R, double *cell_mpCoeff, double *FFCoeff, 
		  double *E, int *Ktable, int2 dof, int2 cutoff, double *VT, 
		  double *Kweights, int grid_type) {

  /* TODO:
   * The current code computes 'CompCoeff' each times for the
   * same 'cell_mpCoeff' when  the cell is in the interaction list of multiple cells. And
   * one optimization is to overwrite 'cell_mpCoeff' with 'CompCoeff' and set the flag,
   * so the next time 'CompCoeff' can be used directly
   */
     
  int n3 = n*n*n, cutoff_f = cutoff.f, cutoff_s = cutoff.s;
  int cutoff2  = cutoff_f * cutoff_s;
  int dofn3    = dof.s * n*n*n;
                    
  // Final multipole expansion (Omega w) =  Sw ; v * Sw = Wl
  int l = 0, l1, l2;
  double tmp, Sw[dofn3];
  double CompCoeff[cutoff_s];

  // 3a in page 8718
  int incr = 1;
  double alpha = 1, beta = 0;
  char trans = 'n';

  // Compute the field proxy values: Matrix Vector multiplication in
  // BLAS: dgemv - 3b
  // cutoff: cutoff on the number of SVD values, 
  // Ecell : M2L operator - E[count] : precalculated, P input, M
  // expansion Output : Pf L expansion

  // Determine the corresponding index in the lookup table
  int k1 = (int)round(R[0]) + 3;
  int k2 = (int)round(R[1]) + 3;
  int k3 = (int)round(R[2]) + 3;

  // note that k1, k2 and k3 are integers from 0 to 6.
  int ninteract = Ktable[49*k1+7*k2+k3];
  assert(ninteract != -1);
  int count; 
     
  if (grid_type == CHEB) {
     
    for (l1=0;l1<n3;l1++) {
      tmp = Kweights[l1];
      for (l2=0;l2<dof.s;l2++) {
	Sw[l] = cell_mpCoeff[l] * tmp;
	l++;
      }
    }

    //print_array(cell_mpCoeff, dof.s*n*n*n, "cell_mpCoeff");
    //print_array(Kweights, n*n*n, "Kweights");
    //print_array(Sw, dof.s*n*n*n, "Sw");
	  
    dgemv_(&trans, &cutoff_s, &dofn3, &alpha, VT, &cutoff_s, Sw, &incr,
	   &beta, CompCoeff, &incr);

    //print_array(CompCoeff, cutoff_s, "CompCoeff");
		  
    count = ninteract*cutoff2;
    dgemv_(&trans,&cutoff_f,&cutoff_s,&alpha,E+count,&cutoff_f,CompCoeff,&incr,
	   &beta,FFCoeff,&incr); // 3b  Ecell is the kernel

    //print_array(FFCoeff, cutoff_f, "FFCoeff");
		  
  } else {
	
    // entry-wise product of cell_mpCoeff
    int N = (int)round(pow(2*n-1, 3));
    int f, s, shift1, shift2=0, shift3=0;
    count = ninteract*(2*n-1)*(2*n-1)*(2*n-1)*dof.s*dof.f;
    for (f=0; f<dof.f; f++, shift2+=N)
      for (s=shift1=0; s<dof.s; s++) {
	FrequencyProduct(N, &E[count+shift3], &cell_mpCoeff[shift1],
			 &FFCoeff[shift2]);
	shift1 += N;
	shift3 += N;
      }	
  }
}


/*
 * Function: Downward Pass
 *  Input:
 *	
 *   	   A: Pointer to the pointer of a box in the hierarchy tree so *A is 
 *            the pointer to the box (see below for variables used in the structure,
 *            refer to file bbfmm.h for the nodeT structure)
 *
 * (*A)->fieldval: Field values at the chebyshev nodes (refer to f on the
 *                 right hand side in step 4 at the bottom of page 8716)
 *
 * (*A)->fieldlist: Pointer to the array of global indices of fields in the box,
 *                  used with 'field' below to find the position of fields in the box
 *
 * (*A)->center: Position of the box center 
 *
 * (*A)->length: Length of the box
 *
 * (*A)->leaves: Pointer to the children boxes, used to determine whether A
 *               is a leaf by checking (*A)->leaves[0] == NULL, and recursivly go down the tree
 *
 * (*A)->leaves[i]->Nf: Number of field points in the child box, used to determine
 *                      whether the it is an empty box to avoid computation by checking
 *                      (*A)->leaves[i]->Nf != 0
 *
 *     field: Pointer to the array of the position of all field points (to be used in L2L) 
 *
 *  Cweights: Pointer to the array of Chebyshev interpolation coefficients, 
 *  	      transfering the weights at Chebyshev nodes between children and parent
 *            (Refer to Sn on the right hand side in step 2 and step 4 
 *            at bottom of page 8715)
 *
 *       Tkz: Pointer to the n x n array (n is the number of chebyshev nodes) 
 *            containing the evaluation of n Chebyshev polynomials 
 *            (T_0 to T_{n-1}) at n chebyshev nodes (of T_n)
 *
 *    cutoff: cutoff.f is the cutoff of the left singular vectors U and 
 *            cutoff.s is the cutoff of the right singular vectors V
 *    
 *	   n: Number of Chebyshev nodes of every dimension in every box 
 *
 *       dof: dof.f is the degree of the field and dof.s is the degree
 *            of the source
 *
 *   Homogen: Order of the homogeneous kernel
 *
 * Output:
 *
 *	phi: Stress at the field points (to be computed in L2P)
 *
 * (*A)->leaves[j]->fieldval: Field values at Chebyshev nodes of children 
 *                            boxes (refer to f on the left hand side in 
 *                            step 4 at the bottom of page 8716)
 *
 * ---------------------------------------------------------------------
 * Distributes the field from the parent cell to the children cells using
 * Chebyshev interpolation. (downward pass of BBFMM)  3c already done. 
 * Only differnece with Upward pass.
 */
void DownwardPass(nodeT **A, vec3 *field, FMMSrc fmm_src, 
		  double *Cweights, double *Tkz,
		  int n, int2 dof, double alpha, kfun_t kfun,
		  double *phi, int grid_type) {

  
  /* Add the contributions from the parent cell to each child cell -
   * otherwise compute all direct interactions and then interpolate 
   * to the field points */
  if ((*A)->leaves[0] != NULL) {

    // Field values for cell 
    double *F = (*A)->fieldval;

    // Determine which children cells contain field points
    int i, l;
    for (i=0;i<8;i++) {
      if ((*A)->leaves[i]->Nf != 0) { 

	double *Fchild = (*A)->leaves[i]->fieldval;

	// Manually set the vector from child to parent
	int  idx[3];
	double r[3];
	idx[2] = i%2, idx[1] = i/2%2, idx[0] = i/4%2;
	for (l=0; l<3; l++)
	  if (idx[l] == 0)
	    r[l] = -1;
	  else
	    r[l] =  1;
		    
	Local2Local(n, r, F, Fchild, dof.f, Cweights, grid_type);
	  
	DownwardPass(&((*A)->leaves[i]), field, fmm_src, Cweights, Tkz,
		     n, dof, alpha, kfun, phi, grid_type);                                                         
      }
    }   
		
  } else { // This is a leaf of the FMM tree.
	
    Local2Particle(A, field, Tkz, n, dof.f, alpha, phi, grid_type);


    /* Due to near field interactions */
    if (fmm_src.src_t == PNT) {

      vec3   *source = fmm_src.source;
      double *q      = fmm_src.q;

  
      int  i, j, k, m, l;
      int  Nf = (*A)->Nf;
      int  nneigh      = (*A)->ineigh;
      int *fieldlist   = (*A)->fieldlist;
      vec3   *fieldpos = malloc(         Nf * sizeof(vec3)   );
      double *fieldval = malloc( dof.f * Nf * sizeof(double) );

      
      // Obtain the positions of the field points
      for (i=0;i<Nf;i++) {
	k = fieldlist[i];
	fieldpos[i].x = field[k].x;
	fieldpos[i].y = field[k].y;
	fieldpos[i].z = field[k].z;
      }

      for (m=0;m<nneigh;m++) {
	nodeT *B        = (*A)->neighbors[m];
	vec3 cshift     = (*A)->cshiftneigh[m];
	int  Ns         = B->Ns;
	int *sourcelist = B->sourcelist;
	vec3 *sourcepos = malloc(       Ns * sizeof(vec3)   );
	double *qsource = malloc( dof.s*Ns * sizeof(double) );
      
	for (j=0;j<Ns;j++) {
	  l = sourcelist[j];
	  sourcepos[j].x = source[l].x + cshift.x;
	  sourcepos[j].y = source[l].y + cshift.y;
	  sourcepos[j].z = source[l].z + cshift.z;
	  for (k=0;k<dof.s;k++)
	    qsource[dof.s*j+k] = q[dof.s*l+k];
	}

	EvaluateField(fieldpos, sourcepos, qsource, Nf, Ns, dof,
		      kfun, fieldval);
      
	for (i=0;i<Nf;i++) {
	  j = dof.f * fieldlist[i];
	  l = dof.f * i;
	  for (k=0;k<dof.f;k++) {
	    phi[j+k] += fieldval[l+k];
	  }
	}
      
	free(sourcepos);
	free(qsource);
      }

      free(fieldval);
      free(fieldpos);
    } // fi: point source

    
  }
}

/* Function: L2L
 * Input:
 *       n: Number of Chebyshev nodes
 *       r: (r[0], r[1], r[2]) is the vector from the child box to the parent box
 *       F: Field values of the parent
 *
 * Output:
 *  Fchild: Field values of the child
 *
 */
void Local2Local(int n, double *r, double *F, double *Fchild, int dof,
		 double *Cweights, int grid_type){

  int l, l1, l2, l3, l4, count1, count2, count3, wstart;
	
  int n2 = n*n;                       // n2 = n^2
  int n3 = n*n*n;                     // n3 = n^3
  int dofn  = dof * n;
  int dofn2 = dof * n2;
  int dofn3 = dof * n3;
	
  //double prefac;
  //if (grid_type) prefac = 2.0/(double)n; // Prefactor for Sn
  //else prefac = 1.;
  //double prefac3 = prefac*prefac*prefac;   // prefac3 = prefac^3
	
  double Fx[dofn3], Fy[dofn3];
	        
  // Recover the index of the child box
  int idx[3] = {0};
  if (  r[0] > 0)
    idx[0] = 1;
  if (  r[1] > 0)
    idx[1] = 1;
  if (  r[2] > 0)
    idx[2] = 1; 

  // Interpolate the parent field along the x-component
  //j = xindex;
  l = 0;
  if (idx[0] == 0)
    wstart = 0;
  else
    wstart = n2;
		
  for (l1=0;l1<n;l1++) {
    count2 = wstart + l1;
    for (l2=0;l2<n;l2++) {
      for (l3=0;l3<n;l3++) {
	count1 = l2*n + l3;
	for (l4=0;l4<dof;l4++) {
	  count3 = dof*count1 + l4;
	  Fx[l]  = ddot_(&n,&F[count3],&dofn2,&Cweights[count2],&n);
	  l++;
	}
      }
    }
  }
		
  // Interpolate the parent field along the y-component
  //j = yindex;
  l = 0;
  if (idx[1] == 0)
    wstart = 0;
  else
    wstart = n2;
			
  for (l1=0;l1<n;l1++) {
    for (l2=0;l2<n;l2++) {
      count2 = wstart + l2;
      for (l3=0;l3<n;l3++) {
	count1 = l1*n2 + l3;
	for (l4=0;l4<dof;l4++) {
	  count3 = dof*count1 + l4;
	  Fy[l] = ddot_(&n,&Fx[count3],&dofn,&Cweights[count2],&n);
	  l++;
	}
      }
    }
  }
		
  /* Interpolate the parent field along the z-component and add
     to child field */
  //j = zindex;
  l = 0;
  if (idx[2] == 0)
    wstart = 0;
  else
    wstart = n2;
			
  for (l1=0;l1<n2;l1++) {
    count1 = l1*n;
    for (l3=0;l3<n;l3++) {
      count2 = wstart + l3;
      for (l4=0;l4<dof;l4++) {
	count3 = dof*count1 + l4;
	Fchild[l] += ddot_(&n,&Fy[count3],&(dof), &Cweights[count2],&n);
	l++;
      }
    }
  }

}

/*Function: L2P
 * Input:
 *
 *     	   A: Pointer to the pointer of a box in the hierarchy tree so *A is 
 *            the pointer to the box (see below for variables used in the structure,
 *            refer to file bbfmm.h for the nodeT structure)
 *
 * (*A)->Nf: Number of fields in the box
 *
 * (*A)->fieldlist: Pointer to the array of global indices of fields in the box,
 *                  used with 'field' below to find the position of fields in the box
 *
 * (*A)->center: Position of the box center 
 *
 * (*A)->length: Length of the box
 *
 * (*A)->fieldval: Filed values at the Chebyshev nodes (refer to the left hand side
 *                 of step 3c at page 8719) 
 *
 *     field: Pointer to the array of the position of all field points 
 *
 *       Tkz: Pointer to the n x n array (n is the number of chebyshev nodes) 
 *            containing the evaluation of n Chebyshev polynomials 
 *            (T_0 to T_{n-1}) at n chebyshev nodes (of T_n)
 *
 *	   n: Number of Chebyshev nodes of every dimension in every box 
 *
 *       dof: dof.f is the degree of the field and dof.s is the degree
 *            of the source
 *
 * Output:
 *
 *       phi: Stress at the field points
 */
 
// Note: no near field interaction
void Local2Particle(nodeT **A, vec3 *field, double *Tkz, int n, int
		    dof, double alpha, double *phi, int grid_type){

  double *F       = (*A)->fieldval;         // Field values for cell 
  vec3 *Sf;
  vec3 *fieldt;                          // Chebyshev-transformed coordinates
  vec3 fcenter = (*A)->center;           // Center of field cell

  int Nf = (*A)->Nf;
  int *fieldlist  = (*A)->fieldlist;        // Indices of field points in cell


  int i, k, l1, l2, l3, idof, lhs_idx, rhs_idx;
  double  tmp, weight;

  double L = AdjustBoxSize((*A)->length, alpha);
  double halfL = L / 2;

	  


  // Map all of the field points to the scaled box
  fieldt   = (vec3 *)malloc(Nf * sizeof(vec3));
  Sf       = (vec3 *)malloc(n * Nf * sizeof(vec3));
  for (i=0;i<Nf;i++) {
    k = fieldlist[i];
    fieldt[i].x = (field[k].x - fcenter.x) / halfL; 
    fieldt[i].y = (field[k].y - fcenter.y) / halfL; 
    fieldt[i].z = (field[k].z - fcenter.z) / halfL; 
  }
		
  // Compute Sf, the mapping function for the field points
  ComputeSn(fieldt, Nf, Sf, n, Tkz, grid_type);
		
  // compute the values due to the far field
  for (i=0; i<Nf; i++) {
    lhs_idx = dof * fieldlist[i];
    rhs_idx = 0;
    for (l1=0; l1<n; l1++) {
      for (l2=0; l2<n; l2++) {
	tmp = Sf[l1*Nf+i].x * Sf[l2*Nf+i].y;
	for (l3=0; l3<n; l3++) {
	  weight = tmp * Sf[l3*Nf+i].z;
		    
	  for (idof=0; idof<dof; idof++)
	    phi[lhs_idx+idof] += weight * F[rhs_idx++];

	}
      }
    }
  }

	    
	  

  /*
  // Map all of the field points to the scaled box
  fieldt   = (vec3 *)malloc(  Nf*numTest * sizeof(vec3));
  Sf       = (vec3 *)malloc(n*Nf*numTest * sizeof(vec3));
  for (i=0;i<Nf*numTest;i++) {
  fieldt[i].x = (fieldRan[i].x - fcenter.x) / halfL; 
  fieldt[i].y = (fieldRan[i].y - fcenter.y) / halfL; 
  fieldt[i].z = (fieldRan[i].z - fcenter.z) / halfL; 
  }
		
  // Compute Sf, the mapping function for the field points
  ComputeSn(fieldt, Nf*numTest, Sf, n, Tkz, grid_type);
		
  // compute the values due to the far field
  for (i=0; i<Nf*numTest; i++) {
  lhs_idx = dof * i;
  rhs_idx = 0;
  for (l1=0; l1<n; l1++) {
  for (l2=0; l2<n; l2++) {
  tmp = Sf[l1*Nf*numTest+i].x * Sf[l2*Nf*numTest+i].y;
  for (l3=0; l3<n; l3++) {
  weight = tmp * Sf[l3*Nf*numTest+i].z;
		    
  for (idof=0; idof<dof; idof++)
  phi[lhs_idx+idof] += weight * F[rhs_idx++];

  }
  }
  }
  }
  */
	  

	  
  free(Sf);
  free(fieldt);
}


/*
 * Function: FreeNode
 * ------------------------------------------------------------------
 * Frees up space associated with node A.
 *
 */
void FreeNode(nodeT *A) {
  int i;
	
  // Free all child nodes first
  for (i=0;i<8;i++) {
    if (A->leaves[i] != NULL) {
      FreeNode(A->leaves[i]);
    }
  }
	
  // Free the arrays for the field and source values
  if (A->fieldval != NULL)
    free(A->fieldval), A->fieldval=NULL;
  if (A->sourceval != NULL)
    free(A->sourceval), A->sourceval=NULL;
  if (A->proxysval != NULL)
    free(A->proxysval), A->proxysval=NULL;
  if (A->sourcefre != NULL)
    fftw_free(A->sourcefre), A->sourcefre=NULL;

  // Free the field and source lists
  if (A->fieldlist != NULL)
    free(A->fieldlist), A->fieldlist=NULL;
  if (A->sourcelist != NULL)
    free(A->sourcelist), A->sourcelist=NULL;
	
  // Last free the node
  free(A);
}

/*
 * Function: ComputeTk
 * ------------------------------------------------------------------
 * Computes T_k(x) for k between 0 and n-1 inclusive.
 */
/*
  void ComputeTk(double x, int n, double *vec) {
  int k;
	
  vec[0] = 1;
  vec[1] = x;
	
  for (k=2;k<n;k++)
  vec[k] = 2.0*x*vec[k-1] - vec[k-2];
  }
*/

/*
 * Function: ComputeSn
 * -----------------------------------------------------------
 * Sn(xm, xi) = -1/n + 2/n* sum_0^(n-1) Tk(xm) Tk(xi) for Chebyshev polynomial
 *
 */
double ClenshawSum(int n, double x, double *Tk);
double LagrangeWeight(int n, double x, int m);
void ComputeSn(vec3 *point, int N, vec3 *Sn, int n, double *Tkz, int grid_type) {
     
  int i, k, m;

  if (grid_type==CHEB) {
	    
    double pfac = 2./n;	    
    double *Tkz_m;

    // Loop over Chebyshev node m
    for (m=0;m<n;m++) {
      k    = m*N;
      Tkz_m = Tkz + m*n;
      for (i=0;i<N;i++) {
	// Compute S_n for each direction using Clenshaw
	Sn[k+i].x = pfac * ClenshawSum(n, point[i].x, Tkz_m);
	Sn[k+i].y = pfac * ClenshawSum(n, point[i].y, Tkz_m);
	Sn[k+i].z = pfac * ClenshawSum(n, point[i].z, Tkz_m);
      }
    }

	    
  } else if (grid_type==UNIF){
	 
    for (m=0;m<n;m++) {
      k = m*N;
      for (i=0;i<N;i++) {
	Sn[k+i].x = LagrangeWeight(n, point[i].x, m);
	Sn[k+i].y = LagrangeWeight(n, point[i].y, m);
	Sn[k+i].z = LagrangeWeight(n, point[i].z, m);
      }
    }
	    
  } // end else
} // end function

	
/*
 * Function: ComputeGrids
 * ---------------------------------------------------------
 * Calculates node locations for grids between [1,-1] and
 * stores them in pre-allocated array nodes.
 */                          
void GridPos1d(double alpha, double len, int n, int grid_type,
	       double* nodes) {
	  
  int m;
  double L = AdjustBoxSize(len, alpha);
  if (grid_type == CHEB) {
    double pi = M_PI;
    for (m=0; m<n; m++)
      nodes[m] = cos(pi*((double)m+0.5)/(double)n) * L;
  }

  else if (grid_type == UNIF)
    for (m=0; m<n; m++)
      nodes[m] = 1 - 2*(double)m/((double)(n-1)) * L;

  else assert(false); // unkown grid type
}

	
// Compute grid points in the box
void GridPos3d(const vec3 center, const double alpha, const double L, const int n, const gridT grid_type, vec3 *P) {

  double gridpos[n]; // [-1, 1]
  //GridPos1d(n, alpha, grid_type, gridpos);
  GridPos1d(alpha, 1.0, n, grid_type, gridpos); 

  int l1, l2, l3, count = 0;
  vec3 vtmp;
  for (l1=0;l1<n;l1++) {
    vtmp.x = (center.x + gridpos[l1]/2) * L;
    for (l2=0;l2<n;l2++) {
      vtmp.y = (center.y + gridpos[l2]/2) * L;
      for (l3=0;l3<n;l3++) {
	P[count].x = vtmp.x;
	P[count].y = vtmp.y;
	P[count].z = (center.z + gridpos[l3]/2) * L;
	count++;
      }
    }
  }

}

	 
/* Get a pair of source and field points with a distance of r units, where r is in [0, 2n-2].
 * Used in ComputeKernelUnif().
 */
void GetPosition(int n, int r, double *fieldpos, double *sourcepos, double *nodepos) {

  if (r < n) {
    *fieldpos  = nodepos[n-1]/2;
    *sourcepos = nodepos[r  ]/2;
  } else {
    *fieldpos  = nodepos[2*(n-1)-r]/2;
    *sourcepos = nodepos[n-1]/2;
  }

}


/*
 * Compute the entry-wise product of two frequencies from rfftw
 * Note 'res' has been initialized
 */
void FrequencyProduct(int N, double *Afre, double *xfre, double *res) {

  int i;
  res[0] += Afre[0]*xfre[0];
  for (i=1; i<N; i++) {
    if (i<(N+1)/2)
      res[i] += Afre[i]*xfre[i] + Afre[N-i]*xfre[N-i];
    else
      res[i] += Afre[N-i]*xfre[i] - Afre[i]*xfre[N-i]; 
  }
  if (N%2 == 0)
    res[N/2] += Afre[N/2]*xfre[N/2];
    
}

/*
// Naive implementation
void Unif2Cheb(double *U, double *C, int n, int dof) {

double Ugrid[n];
GridPos1d(n, 0, UNIF, Ugrid);
//ComputeGrids(n, Ugrid, UNIF);
//ComputeGrids(n, Ugrid, CHEB);

// 3d uniform grids
int l1, l2, l3, count, n3=n*n*n;
vec3 Ugrid_3d[n3];
count = 0;
for (l1=0; l1<n; l1++)
for (l2=0; l2<n; l2++)
for (l3=0; l3<n; l3++) {
Ugrid_3d[count].x = Ugrid[l1];
Ugrid_3d[count].y = Ugrid[l2];
Ugrid_3d[count].z = Ugrid[l3];
count++;
}
  
vec3 U2C[n3*n];
ComputeSn(Ugrid_3d, n3, U2C, n, Tkz, CHEB);
//ComputeTk(Tkz, n);
//ComputeSn(Ugrid_3d, Tkz, n, n3, U2C, UNIF);

int l4, j, k;
double sum;
count = 0;
for (l1=0; l1<n; l1++) 
for (l2=0; l2<n; l2++) 
for (l3=0; l3<n; l3++)
for (l4=0; l4<dof; l4++) {
sum = 0;
for (j=0; j<n3; j++) {
k = dof * j + l4;
sum += U[k]*
U2C[l1*n3+j].x*
U2C[l2*n3+j].y*
U2C[l3*n3+j].z;
}
C[count++] = sum;
}
}


// Naive algorithem
void Cheb2Unif(double *C, double *U, int n, int dof) {

double Cgrid[n];
GridPos1d(n, 0, UNIF, Cgrid);
//ComputeGrids(n, Cgrid, CHEB);
//ComputeGrids(n, Cgrid, UNIF);

// 3d Chebyshev grids
int l1, l2, l3, count, n3=n*n*n;
vec3 Cgrid_3d[n3];
count = 0;
for (l1=0; l1<n; l1++)
for (l2=0; l2<n; l2++)
for (l3=0; l3<n; l3++) {
Cgrid_3d[count].x = Cgrid[l1];
Cgrid_3d[count].y = Cgrid[l2];
Cgrid_3d[count].z = Cgrid[l3];
count++;
}

  
vec3 C2U[n3*n];
//ComputeSn(Cgrid_3d, n3, n, C2U, UNIF);
ComputeSn(Cgrid_3d, n3, C2U, n, Tkz, CHEB);


int i, k, l4, l;
double sum = 0, tmp1, tmp2;
for (i=0; i<n3; i++) {
k = dof * i;
for (l4=0; l4<dof; l4++) {
sum = 0;
l = l4;
for (l1=0; l1<n; l1++) {
tmp1 = C2U[l1*n3+i].x;
for (l2=0; l2<n; l2++) {
tmp2 = tmp1 * C2U[l2*n3+i].y;
for (l3=0; l3<n; l3++) {
sum += C[l] * tmp2 * C2U[l3*n3+i].z;
l   += dof;
}
}
}
U[k++] = sum;
}
}
  
}
*/
 
// Evaluate Chebyshev polynomial form order 0 to n-1 at n Chebyshev nodes. Tk: nodes x order  in 1d array
void ComputeTkz(double *Tkz, int n) {
   
  // Compute n grids in [-1, 1]
  double nodes[n], vec[n];
  GridPos1d(0, 1.0, n, CHEB, nodes);
	
  int i, k, m;
  double x;
  for (m=0; m<n; m++) {

    x = nodes[m];
     
    vec[0] = 1;
    vec[1] = x;
	
    for (k=2;k<n;k++)
      vec[k] = 2.0*x*vec[k-1] - vec[k-2];

    i = m*n;
    for (k=0;k<n;k++)
      Tkz[i+k] = vec[k];
  }
}


// summation using Clenshaw's recurrence relation
double ClenshawSum(int n, double x, double *Tk) {
  int j;
  double d0, d1, d2;
  d0 = d1 = 0;
  for (j = n - 1; j > 0; j--) {
    d2 = d0;
    d0 = 2.0 * x * d0 - d1 + Tk[j];
    d1 = d2;
  }
  return x * d0 - d1 + 0.5 * Tk[0];
}


double LagrangeWeight(int n, double x, int m) {
  int j;
  double num = 1., denom = 1.;
  for (j=0;j<n;j++) {
    if(m!=j) {
      num   *= x - (1-2*(double)j/(double)(n-1));
      denom *= -2*((double)m-(double)j)/(double)(n-1);
    }
  }
  return num/denom;
}
 
 
double AdjustBoxSize(double L, double alpha){
  return L * (1+alpha);
}


// check if two boxes are well seperated
bool Is_well_separated(vec3 p, double L) {
  /*
  int x = round(abs((1.0 / L) * p.x));
  int y = round(abs((1.0 / L) * p.y));
  int z = round(abs((1.0 / L) * p.z));
*/
    
  int x = round(fabs(p.x / L));
  int y = round(fabs(p.y / L));
  int z = round(fabs(p.z / L));
  return x > 1 || y > 1 || z > 1;
}
 

void Create_lookup_table(int *Ktable) {
  int l1, l2, l3;
  int ncell = 0, ninteract = 0;
  for (l1=-3;l1<4;l1++) {
    for (l2=-3;l2<4;l2++) {
      for (l3=-3;l3<4;l3++) {
	vec3 r={l1, l2, l3};
	if ( Is_well_separated(r, 1.0) ) {
	  Ktable[ncell] = ninteract;
	  ninteract++;
	}
	else Ktable[ncell] = -1;
	    
	ncell++;
      }
    }
  }	
}



// Note: Out should be initialized first
void Interpolate(double *In, int in_gridtype, double *Out, int out_gridtype, int n, int dof, double *Tkz) {

  int n3 = n*n*n;
  vec3 out_grid[n3], center={0,0,0};
  GridPos3d(center, 0, 2, n, out_gridtype, out_grid);

  vec3 Sn[n3*n];
  ComputeSn(out_grid, n3, Sn, n, Tkz, in_gridtype);

  int i, l1, l2, l3, idof;
  int lhs_idx, rhs_idx;
  double w1, w2, w3;

  for (i=0; i<n3; i++) {

    lhs_idx = dof * i;
    rhs_idx = 0;
    for (l1=0; l1<n; l1++) {
      w1 = Sn[l1*n3+i].x;
      for (l2=0; l2<n; l2++) {
	w2 = w1 * Sn[l2*n3+i].y;
	for (l3=0; l3<n; l3++) {
	  w3 = w2 * Sn[l3*n3+i].z;
		    
	  for (idof=0; idof<dof; idof++)
	    Out[lhs_idx+idof] += w3 * In[rhs_idx++];

	}
      }
    }

  }
}


// Note: Out should be initialized first
void Anterpolate(double *In, int in_gridtype, double *Out, int out_gridtype, int n, int dof, double *Tkz) {

  int n3 = n*n*n;
  vec3 in_grid[n3], center = {0,0,0};
  GridPos3d(center, 0, 2, n, in_gridtype, in_grid);
  
  vec3 Sn[n3*n];
  ComputeSn(in_grid, n3, Sn, n, Tkz, out_gridtype);


  int l1, l2, l3, l4, j;
  int x_idx, y_idx, z_idx, lhs_idx, rhs_idx;
  double weight;

  lhs_idx = 0;
  for (l1=0;l1<n;l1++) {
    x_idx = l1*n3;
    for (l2=0;l2<n;l2++) {
      y_idx = l2*n3;
      for (l3=0;l3<n;l3++) {
	z_idx = l3*n3;

	for (j=0; j<n3; j++) {
	  weight = Sn[x_idx+j].x *Sn[y_idx+j].y *Sn[z_idx+j].z;
	  rhs_idx = j * dof;
		  
	  for (l4=0; l4<dof; l4++)
	    Out[lhs_idx+l4] += weight * In[rhs_idx++];

	}
	lhs_idx += dof;
		
      }
    }
  }
}


double* MeanStressFmm(const double *ChebyshevWeightField, const int n, const int dof, const double *Tkz) {

  // Integrate Sn
  int i, j, l1, l2, l3;
  double SnInt[n];
  for (i=0; i<n; i++) {
    SnInt[i] = 1;
    for (j=2; j<n; j+=2)
      // Int_{-1}^1 Tk(x) dx
      SnInt[i] -= 2./(j*j-1) * Tkz[i*n+j]; 
    
    SnInt[i] *= 2./n;
  }  
      
  // Compute mean stress
  int idof, count;
  double meanWeight;
  double *meanStressFmm = calloc( dof, sizeof(double) );      
  for (l1=0; l1<n; l1++)
    for (l2=0; l2<n; l2++)
      for (l3=0; l3<n; l3++) {
	count = (l1 * n*n + l2*n + l3) * dof;
        
	meanWeight = SnInt[l1] * SnInt[l2] * SnInt[l3];
        
	for (idof=0; idof<dof; idof++)    
	  meanStressFmm[idof] += ChebyshevWeightField[count + idof] * meanWeight;
      }
  
  // divide by 8 that is the volume
  meanStressFmm = Multiply (meanStressFmm, 1./8, dof);
  return meanStressFmm;
}

 
bool IsHomoKernel( double homogen ) {
  return fabs(homogen) > HOMO_THRESHOLD;
}


// Compute the midpoint of every segment
vec3* ComputeMiddlePoint( segT *segment, int Ns ) {

  vec3 *midpoint = (vec3 *) malloc(Ns * sizeof(vec3));
  
  int i;
  for (i=0; i<Ns; i++) {
    midpoint[i].x = (segment[i].p_beg.x + segment[i].p_end.x)/2;
    midpoint[i].y = (segment[i].p_beg.y + segment[i].p_end.y)/2;
    midpoint[i].z = (segment[i].p_beg.z + segment[i].p_end.z)/2;
  }

  return midpoint;
}



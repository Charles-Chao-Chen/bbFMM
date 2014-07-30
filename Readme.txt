TODO: consider more symmetry in Chebyshev pre-computation.
Only the r- symmetry is currently considered, 
i.e. f(r) ?= f(-r).

Changed 'FMMSetup()' to 'ComputeM2L()'.

Put 'FMMReadMatrices()' into 'FMMSetup()'.

Created 'BuildTree()'.

Changed 'dof2' to 'int2'.

Replace fftw2 with fftw3. Done 

Set uniform 'Umat' and 'Vmat' to be existing files
in FMMSetup().

uniform grids for non-homogeneous kernel. Done

Removed skipping level in non-homogeneous kernel

Todo: use macro for read number check. Done

Bug: Gaussian set for 54, but dof={1,1} in main.c

Adjust box size in PBC to be tested.

ComputeWeightsPBC() to be optimized.

U = VT = NULL in Uniform grids

Gave up merging Interpolate() and Local2Particle();
Anterpolate() and Particle2Moment(), because of the
annoying 'sourcelist' and 'fieldlist'.

Moved DirPBC() and DirWo1() from unility to directCalc;
Fixed memory lost issue in AddPBC().

Use Interpolate() and Anterpolate()

TODO: implement U(l1,l2,l3) = sum S(l1,k1) S(l2,k3) S(l3,k3) C(k1,k2,k3)

TODO: implement boxLen adjustment

utility.h

Return nodes in [-.5, .5]

Use symmetry in ComputeKernelUniform() to save a half.

Move 'prefac = 2./n' into 'ComputeSn()'

Separate kernel definition in one file.

PBC problem discovered in Parameters->e3.

Discover 'Ksize' in bbfmm() overflow problem.

Uniform grids is tested, but the PBC gives ignorant result.

In P2M, dof loop should be inside. L2P change loop order
to have contiguous memory access.

Optimize multi-rhs

Fix DirctCalc3D() for multi-dof.

This version modifies for the non-singular when source=field
in EvaluateField() by detecting infinity value.
e.g. kernel = exp(-x). 

This version adds unifrom grids with fft. 
TODO: PBC, adjusting box size, non-symmetric,
non-homogeneous kernels.

This version implement non-symmetric parameter 'symmtry'
and non-homogeneous with the existing parameter homogen=0.
TODO: Implement the case where new pre-computation
contains the exsiting data.
e.g.: fileBoxLen = 1, fileTreeLevel = 2
and boxLen = 2, fileTreeLevel = 3.

This version tests PBC results.

This version use a faster way to compute aniso-
tropic kernel and test accuracy (0 + pbc)

This version has a clean interface for computing
stress within the original domain plus PBC stress 
(subtract mean stress).

Note: ComputeKernelSVD() only works for anti-symmetry
kernels

*Note: only LINEINT version is available

This version calculates mean stress from PBC.

This version adjust the box size to accomodate 
segments outside boxes.

This version creates M2M,  L2L and M2L routines needed 
in ParaDis.

The M2M function takes in the source value of one 
children and returns its contribution to the parent.

The L2L function takes in the field value of the parent
and returns the anterpolatino value of the children.

The M2L function computes the contribution from one
cell's multiple coefficients and get the far field
interaction. 

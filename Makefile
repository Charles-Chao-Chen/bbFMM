CC = gcc
LD = gcc

#BLAS_DIR = /cm/shared/apps/openblas/dynamic/0.2.6/lib
#BLAS_INCLUDE = /cm/shared/apps/openblas/dynamic/0.2.6/include
#LAPACK_DIR = /cm/shared/apps/lapack/open64/64/3.4.2
#FFTW_DIR = /cm/shared/apps/fftw/openmpi/open64/64/2.1.5/double/lib
#FFTW_INCLUDE = /cm/shared/apps/fftw/openmpi/open64/64/2.1.5/double/include

# make sure the corresponding modules have been added
# and use 'extern "C"{}' to declare lapack routines
LDPATH  = -L /usr/bin/ld
LDFLAGS = -lblas -llapack -lfftw3 -lm 
FLAGS   = -Wall -O3 -g -I $(BLAS_INCLUDE) -I $(FFTW_INCLUDE)
#KERNEL  = -DDDANISO -DLINEINT
#KERNEL = -DLAPLACIANFORCE -DTENSOR
#KERNEL_INT   = -DLAPLACIANFORCE -DLINEINT
KERNEL  = -DTENSOR #-DLINEINT
PFLAG   =


PROG   = bbfmm
HEADER = AnisoFunctions.h bbfmm.h kernelFun.h directCalc.h utility.h
OBJS   = main.o bbfmm.o directCalc.o utility.o kernelFun.o DerivativeOfGreensFunction.o ChaosFunction.o AnisoFunctions.o

all: $(PROG)

bbfmm: $(OBJS)
	$(LD) $(OBJS) -o $@ $(LDPATH) $(LDFLAGS) $(PFLAG)


main.o: main.c $(HEADER)
	$(CC) $(FLAGS) $(KERNEL) -c $< $(PFLAG)


bbfmm.o: bbfmm.c $(HEADER)
	$(CC) $(FLAGS) $(KERNEL) -c $< -o $@ $(PFLAG)


kernelFun.o: kernelFun.c $(HEADER)
	$(CC) $(FLAGS) $(KERNEL) -c $< -o $@ $(PFLAG)

directCalc.o: directCalc.c $(HEADER)
	$(CC) $(FLAGS) $(KERNEL) -c $< -o $@ $(PFLAG)

utility.o: utility.c $(HEADER)
	$(CC) $(FLAGS) $(KERNEL) -c $< -o $@ $(PFLAG)

DerivativeOfGreensFunction.o: DerivativeOfGreensFunction.c
	$(CC) $(FLAGS) -Wno-unused-variable -c $< $(PFLAG)

ChaosFunction.o: ChaosFunction.c
	$(CC) $(FLAGS) -c $< $(PFLAG)

AnisoFunctions.o: AnisoFunctions.c AnisoFunctions.h
	$(CC) $(FLAGS) -c $< $(PFLAG)

#%.o:	%.c %.h
#	$(CC) $(FLAGS) -c $< $(PFLAG)

clean:
	/bin/rm -rf *.o *~ \#*\# $(PROG) *core

tar:
	tar cvfz 3DbbFMM.tgz Makefile main.c bbfmm.c bbfmm.h kernelFun.c kernelFun.h directCalc.c directCalc.h utility.c utility.h AnisoFunctions.c AnisoFunctions.h ChaosFunction.c DerivativeOfGreensFunction.c Readme.txt

tags:
	etags *.h *.c

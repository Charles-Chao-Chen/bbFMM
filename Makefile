CC = gcc
LD = gcc

FFTW_ROOT = /home/cchen10/Softwares/fftw-3.3.5
FFTW_LIB = $(FFTW_ROOT)/lib
FFTW_INCLUDE = $(FFTW_ROOT)/include

# make sure the corresponding modules have been added
# and use 'extern "C"{}' to declare lapack routines
LDPATH  = -L /usr/bin/ld
LDFLAGS = -lblas -llapack -L$(FFTW_LIB) -lfftw3 -lm 
FLAGS   = -Wall -O3 -I$(FFTW_INCLUDE)
PFLAG   =


PROG   = bbfmm test
HEADER = AnisoFunctions.h bbfmm.h kernelFun.h directCalc.h utility.h
OBJS   = bbfmm.o directCalc.o utility.o kernelFun.o DerivativeOfGreensFunction.o ChaosFunction.o AnisoFunctions.o

OBJS1  = main.o $(OBJS)
OBJS2  = test.o $(OBJS)


all: $(PROG)

bbfmm: $(OBJS1)
	$(LD) $(OBJS1) -o $@ $(LDPATH) $(LDFLAGS) $(PFLAG)

test: $(OBJS2)
	$(LD) $(OBJS2) -o $@ $(LDPATH) $(LDFLAGS) $(PFLAG)

test.o: test.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)

main.o: main.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)


bbfmm.o: bbfmm.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)


kernelFun.o: kernelFun.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)

directCalc.o: directCalc.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)

utility.o: utility.c $(HEADER)
	$(CC) $(FLAGS) -c $< -o $@ $(PFLAG)

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

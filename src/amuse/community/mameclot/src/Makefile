CC = cc
CFLAGS = -std=c99 -Wall -W -O3
EXE = mameclot
SOURCE = $(EXE).c 
OBJ = $(EXE).o
CPUOBJ = pot.o
LIBS = -lm

CUDAPATH = /usr/local/cuda/bin/
CUDAFLAGS = -L/usr/local/cuda/lib64 -I/usr/local/cuda/include -lcudart -lcuda
NVCC = $(CUDAPATH)/nvcc
CUDASOURCE = pot.cu
CUDAOBJ = gpupot.o


default:	$(OBJ) $(CPUOBJ) 
		$(CC)  -o $(EXE) pot.o $(OBJ) $(LIBS) 

gpu:	$(OBJ)
	$(NVCC) -c -m64 -arch sm_20  gpupot.cu
	$(CC)  $(CUDAFLAGS)  -o $(EXE).gpu $(CUDAOBJ) $(OBJ) $(LIBS)


clean:	
	rm -f *.o $(EXE) $(EXE).gpu



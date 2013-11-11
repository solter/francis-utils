# Francis algorithm in C; see README

CC = gcc
CFLAGS = -I/usr/local/include
CFLAGS += -Wall#-pedantic -Wall
LDFLAGS = -L/usr/lib -L/usr/lib64/atlas
LDLIBS = -llapack -lblas -lm
#-llapack_atlas -latlas
SOURCES = testbench.c bulge.c util.c
ODIR = objDir

test: testbench.x
	./testbench.x

testbench.x: $(patsubst %.c,$(ODIR)/%.o,$(SOURCES)) | $(ODIR)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(ODIR)/%.o: %.c | $(ODIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(ODIR):
	mkdir $@

clean: objDir
	rm -f testbench.x $(patsubst %.c,$(ODIR)/%.o,$(SOURCES)) 

.PHONY: clean .CLEAN

#C=nvcc
#  2 GCC=g++
#  3 #NVCC_OPTS=-O3 -arch=sm_20 -Xcompiler -Wall -Xcompiler -Wextra -m64
#  4 NVCC_OPTS=
#  5 CUSP_LOC=-I cusplibrary/
#  6 EX_LIBS=-lcurand
#  7 CODE=code/
#  8 
#  9 main.x: objDir/main.o objDir/fileOps.o objDir/helper.o objDir/kernels.o
# 10   $(NVCC) -o $@ $^ $(NVCC_OPTS) $(CUSP_LOC) $(EX_LIBS)
# 11 
# 12 objDir/kernels.o: $(CODE)kernels.cu | objDir
# 13   $(NVCC) -c $< -o $@ $(NVCC_OPTS) $(CUSP_LOC)
# 14 
# 15 objDir/main.o: $(CODE)main.cu | objDir
# 16   $(NVCC) -c $< -o $@ $(NVCC_OPTS) $(CUSP_LOC)
# 17 
# 18 objDir/fileOps.o: $(CODE)fileOps.cu | objDir
# 19   $(NVCC) -c $< -o $@ $(NVCC_OPTS) $(CUSP_LOC)
# 20 
# 21 objDir/helper.o: $(CODE)helper.cu | objDir
# 22   $(NVCC) -c $< -o $@ $(NVCC_OPTS) $(CUSP_LOC)
# 23 
# 24 objDir:
# 25   mkdir objDir
# 26 
# 27 test: main.x
# 28   ./main.x 12A.mtx 12I.mtx singTest outs.txt 
# 29 
# 30 clean: objDir
# 31   rm -f objDir/*.o *.bin


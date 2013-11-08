
# Francis algorithm in C; see README

CC = clang
CFLAGS = -I/usr/local/include
CFLAGS += -pedantic -Wall
LDFLAGS = -L/usr/lib -L/usr/lib64/atlas
LDLIBS = -llapack -lblas -lm
#-llapack_atlas -latlas
SOURCES = testbench.c bulge.c util.c

test: testbench
	./testbench

testbench: ${SOURCES:.c=.o}

clean:
	rm -f testbench ${SOURCES:.c=.o}

.PHONY: clean .CLEAN


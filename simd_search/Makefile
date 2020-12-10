CC ?= gcc
ifndef CFLAGS
	CFLAGS := ${CFLAGS} -g -std=c99 -march=native -mtune=native -O3

	# Flags from "Parallel and High Performance Computing" by Robert Robey and Yuliana Zamora
	ifeq ($(CC),gcc)
		CFLAGS := ${CFLAGS} -fopenmp-simd -fstrict-aliasing -ftree-vectorize -mprefer-vector-width=512 -fno-trapping-math -fno-math-errno -fopt-info-vec-all
	else ifeq ($(CC),clang)
		CFLAGS := ${CFLAGS} -fopenmp-simd -fstrict-aliasing -fvectorize -fno-math-errno -Rpass-analysis=loop-vectorize
	else ifeq ($(CC),icc)
		CFLAGS := ${CFLAGS} -qopenmp-simd -ansi-alias -restrict -fp-model:precise -xHOST -vecabi=cmdtarget -qopt-zmm-usage=high -qopt-report=5
	endif
endif

all: test

algorithms.o: algorithms.c algorithms.h
	CC $(CFLAGS) -c algorithms.c

test: test.c algorithms.o
	CC $(CFLAGS) algorithms.o test.c -o test

.PHONY: clean
clean:
	rm algorithms.o test
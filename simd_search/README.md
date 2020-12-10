# SIMD-Optimized Search Algorithms
This directory contains C code for the search algorithms discussed in *SIMD-Optimized Search Over Sorted Data*.

[*Insert link to paper here when published*]

## Usage
Run `make` to build, then `./test` to run. The Makefile automatically chooses the best vectorization flags depending on the compiler you pass. Best vectorization results are found with the Intel compiler: `make CC=icc`. To replace the flags with your own, simply specify `CFLAGS`. For example, if you want to run with no optimization and debug symbols, run `make CFLAGS='-g -O0'`.

## Algorithms
### Logarithm Hash Search
### Skiplist Search
### Exponent Hash Search
### Hunt & Locate Search

## Authors
 - *Benjamin Mastripolito* - bmastripolito@lanl.gov
 - *Nicholas Koskelo* - nkoskelo@lanl.gov
 - *Dylan Weatherred* - dweatherred@lanl.gov
 - *David A. Pimentel* - davidp@lanl.gov
 - *Daniel Sheppard* - danielsheppard@lanl.gov
 - *Anna Pietarila Graham* - annap@lanl.gov
 - *Laura Monroe* - lmonroe@lanl.gov
 - *Robert Robey* - brobey@lanl.gov

## Issues
Sometimes when using all performance flags, search will segfault.
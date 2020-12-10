/**
 * @file test.c
 * @author Benjamin Mastripolito (bmastripolito@lanl.gov)
 * @brief Testing search algorithms from algorithms.c
 * @date 2020-12-09
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "algorithms.h"

// Change to #define to check for validity of searches
#undef CHECK_VALIDITY

typedef struct {
    const char* name;
    int* (*search_fn)(double* x, int n, double* y, int m);
} search_fn_t;

double randreal(double min, double max) {
    return min + ((double)rand() / (double)RAND_MAX) * (max - min);
}

int main(int argc, char const *argv[]) {
    // Search functions to test
    const search_fn_t searchFuncs[] = {
        {"Binary Search",           binarySearch},
        {"Linear Search",           linearSearch},
        {"Hunt & Locate Search",    huntLocateSearch},
        {"Logarithm Hash Search",   logHashSearch},
        {"Skiplist Search",         skiplistSearch},
        {"Exponent Hash Search",    expHashSearch},
    };
    const int nSearchFuncs = sizeof(searchFuncs) / sizeof(search_fn_t);
    // Number of data array elements
    const int n = 100;
    // Number of targets
    const int m = 10000000;
    // The range of powers of 10 between smallest and largest element
    const double expRange = 50.0;
    double* x = (double*)malloc(n * sizeof(double));
    double* y = (double*)malloc(m * sizeof(double));

    // Generate random data
    for (int i = 0; i < n; i++) {
        x[i] = pow(10.0, (((double)i / (double)n) - 0.5) * expRange);
    }
    for (int i = 0; i < m; i++) {
        y[i] = randreal(0.0, 1.0) * pow(10.0, randreal(-expRange * 0.65, expRange * 0.65));
    }

    // Perform a linear search for validation
    int* ref = linearSearch(x, n, y, m);

    // Time search functions
    for (int i = 0; i < nSearchFuncs; i++) {
        struct timeval tStart, tEnd;
        gettimeofday(&tStart, NULL);
        int* compare = searchFuncs[i].search_fn(x, n, y, m);
        gettimeofday(&tEnd, NULL);
        unsigned long usec = ((1000000 * (tEnd.tv_sec - tStart.tv_sec)) + tEnd.tv_usec) - tStart.tv_usec;
        printf("%s took %.5lf seconds\n", searchFuncs[i].name, usec / 1000000.0);
        // Validate results
        #ifdef CHECK_VALIDITY
        for (int j = 0; j < m; j++) {
            if (compare[j] != ref[j]) {
                printf("INVALID %d> target=%le ref=%d comp=%d\n", j, y[j], ref[j], compare[j]);
            }
        }
        #endif
        free(compare);
    }

    free(ref);
    free(x);
    free(y);
    return 0;
}

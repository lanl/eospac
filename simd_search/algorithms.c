/**
 * @file algorithms.c
 * @author Benjamin Mastripolito (bmastripolito@lanl.gov)
 * @brief Definitions for various SIMD-enabled search algorithms
 * @date 2020-12-09
 */

#include "algorithms.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Branchless minimum of two integers
 * 
 * @param a First integer
 * @param b Second integer
 * @return int Smaller of a and b
 */
static inline int min(int a, int b) {
    const int cond = a < b;
    return (cond * a) | (!cond * b);
}

/**
 * @brief Branchless maximum of two integers
 * 
 * @param a First integer
 * @param b Second integer
 * @return int Larger of a and b
 */
static inline int max(int a, int b) {
    const int cond = a > b;
    return (cond * a) | (!cond * b);
}

/**
 * @brief Branchless choice of two integers
 * 
 * @param condition 1 for True, 0 for False
 * @param valTrue Return if condition == 1
 * @param valFalse Return if condition == 0
 * @return int Either valTrue or valFalse
 */
static inline int choose(int condition, int valTrue, int valFalse) {
    return (condition * valTrue) | (!condition * valFalse);
}

/**
 * @brief Inlined binary search of x for lower bound of target between specified indices
 * 
 * @param start Starting index (inclusive)
 * @param end Ending index (inclusive)
 * @param x Array to search within
 * @param target The target value to search for
 * @return int The index of the last element in x that is smaller than target
 */
static inline int inline_binarySearch(int start, int end, double *restrict x, double target) {
    while ((end - start) > 1) {
        const int midPoint = (end + start) / 2;
        const int c = target < x[midPoint];
        end = choose(c, midPoint, end);
        start = choose(!c, midPoint, start);
    }
    return choose(target < x[end], start, end);
}

/**
 * @brief Inlined linear search of x for lower bound of target between specified indices
 * 
 * @param start Starting index (inclusive)
 * @param end Ending index (inclusive)
 * @param x Array to search within
 * @param target The target value to search for
 * @return int The index of the last element in x that is smaller than target
 */
static inline int inline_linearSearch(int start, int end, double *restrict x, double target) {
    int i = start + 1;
    for (; i <= end; i++) {
        if (x[i] > target) break;
    }
    return i - 1;
}

/**
 * @brief Returns log-based hash of value given adjustment parameters
 *
 * @param x The value to calculate the hash of
 * @param baseAdjust 1/log10(B) where B is the base to convert to
 * @param baseOffset Integer value to offset the hash with (it is added to the result)
 * @return The hash result
 */
static inline long logHash(double x, double baseAdjust, long baseOffset) {
    return (long)(log10(x) * baseAdjust) + baseOffset;
}

int *logHashSearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    // Account for negative numbers, adds magic number to avoid zeroes
    double minVal = x[0];
    double offset = 0.0;
    if (minVal < 0.0) {
        offset = -(minVal * 0.012323745192304);
    } else if (minVal == 0.0) {
        offset = x[1] * 0.012323745192304;
    }
    minVal = x[0] + offset;

    // Determine the base we need to convert to, and the baseOffset
    const int htableSize = 512;
    const double maxVal = x[n - 1] + offset;
    const double base = pow(10.0, (log10(maxVal) - log10(minVal)) / (htableSize - 1));
    const double baseAdjust = 1.0 / log10(base);
    const long baseOffset = -(log10(minVal) / log10(base));

    // Populate the hash table
    int hashTable[htableSize];
    #pragma omp simd aligned(hashTable)
    for (int i = 0; i < htableSize; i++) {
        hashTable[i] = n - 1;
    }
    for (int i = 0; i < n; i++) {
        const int index = logHash(x[i] + offset, baseAdjust, baseOffset);
        hashTable[index] = min(hashTable[index], i);
    }

    // Correct for if an order of magnitude got skipped
    for (int k = htableSize - 1; k > 0; k--) {
        const int upperBound = hashTable[k];
        while (hashTable[k - 1] > upperBound) {
            hashTable[k - 1] = upperBound - 1;
            k--;
        }
    }

    // Perform the search
    #pragma omp simd aligned(lowbounds, y)
    for (int i = 0; i < m; i++) {
        const double target = y[i];
        // Uses branching for bounds checking, see expHashSearch for bounds checking without branching
        if (target + offset < minVal) {
            lowbounds[i] = 0;
        } else if (target + offset > maxVal) {
            lowbounds[i] = n - 1;
        } else {
            const long hash = logHash(target + offset, baseAdjust, baseOffset);
            // Retrieve the indices to search between
            const int start = max(hashTable[hash] - 1, 0);
            const int end = choose(hash < htableSize - 1, hashTable[hash + 1], n - 1);
            // Perform a linear search for target, replace with inline_binarySearch if desired
            lowbounds[i] = inline_linearSearch(start, end, x, target);
        }
    }

    return lowbounds;
}

int *skiplistSearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    // Elements per cacheline
    const int epc = 64 / sizeof(double);

    // Skiplist setup
    const int skiplistSize = n / epc;
    double skiplist[skiplistSize];
    for (int i = 0; i < skiplistSize; i++) {
        skiplist[i] = x[min(n - 1, i * epc)];
    }
    const int maxListIndex = (skiplistSize - 1) * epc;
    const double maxListVal = skiplist[skiplistSize - 1];
    const double maxVal = x[n - 1];
    const double minVal = x[0];

    // Perform the search
    #pragma omp simd aligned(lowbounds, y)
    for (int i = 0; i < m; i++) {
        const double target = y[i];
        if (target < minVal) {
            lowbounds[i] = 0;
        } else if (target > maxVal) {
            lowbounds[i] = n - 1;
        } else {
            int start = maxListIndex;
            int end = n - 1;
            if (target < maxListVal) {
                start = inline_binarySearch(0, skiplistSize - 1, skiplist, target) * epc;
                end = min(start + epc - 1, n - 1);
            }
            lowbounds[i] = inline_linearSearch(start, end, x, target);
        }
    }

    return lowbounds;
}

/**
 * @brief Returns the base-2 exponent of x by extracting from the exponent bits
 * 
 * @param x The value to extract the exponent from
 * @return int The exponent
 */
static inline int expHash(double x) {
    // Converts x to an unsigned long pointer so that the exponent (bits 62-52) can be extracted
    // Subtracts 1022 because that is the offset to normalize the stored value to a signed value
    // https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    return (int)(0x7ff & ((*(unsigned long *)(&x)) >> 52)) - 1022;
}

int *expHashSearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    // Account for negative numbers
    double offset = 0.0;
    if (x[0] < 0.0) {
        // Add a number of smaller magnitude so that we don't have to deal with zeroes
        offset = -(x[0] * 0.012323745192304);
    } else if (x[0] == 0.0) {
        offset = x[1] * 0.012323745192304;
    }

    const double minVal = x[0] + offset;
    const double maxVal = x[n - 1] + offset;
    const int minHash = expHash(minVal);
    const int htableSize = expHash(maxVal) - minHash + 1;
    int hashTable[htableSize];

    // Initialize the hash table with highest index into x
    #pragma omp simd aligned(hashTable)
    for (int i = 0; i < htableSize; i++) {
        hashTable[i] = n - 1;
    }

    // Populate the hash table with indices
    #pragma omp simd aligned(x)
    for (int i = 0; i < n; i++) {
        const int index = expHash(x[i] + offset) - minHash;
        hashTable[index] = min(hashTable[index], i);
    }

    // Correct for if an order of magnitude got skipped.
    for (int k = htableSize - 1; k > 0; k--) {
        int upperBound = hashTable[k];
        while (hashTable[k - 1] > upperBound) {
            hashTable[k - 1] = upperBound - 1;
            k--;
        }
    }

    // Perform search
    #pragma omp simd aligned(y, lowbounds)
    for (int i = 0; i < m; i++) {
        const double target = y[i];
        // Branchless clamp the hash value to the bounds of the hash table
        const int hash = expHash(target + offset) - minHash;
        // Retrieve the indices to search between
        const int bounded = hash >= 0 && hash < htableSize;
        const int start = choose(bounded, min(hashTable[hash] - 1, 0), 1);
        const int end = choose(bounded, choose(hash < htableSize - 1, hashTable[hash + 1], n - 1), 0);

        int jlow = inline_linearSearch(start, end, x, target);
        jlow = choose(target < minVal, 0, jlow);
        jlow = choose(target > maxVal, n - 1, jlow);
        lowbounds[i] = jlow;
    }
    return lowbounds;
}

int *huntLocateSearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    const double minVal = x[0];
    const double maxVal = x[n - 1];
    int start = 0;
    int end = n - 1;

    for (int i = 0; i < m; i++) {
        const double target = y[i];
        if (target < minVal) {
            lowbounds[i] = 0;
        } else if (target > maxVal) {
            lowbounds[i] = n - 1;
        } else {
            // "Hunt" for bounds, re-using last bounds
            int j = 1;
            while (start >= 0 && end < n) {
                if (target > x[end]) {
                    start = end;
                    end += j;
                } else if (target < x[start]) {
                    end = start;
                    start -= j;
                } else {
                    break;
                }
                j <<= 1;
            }
            // Binary search within bounds
            end = min(end, n - 1);
            start = inline_binarySearch(max(start, 0), end, x, target);
            lowbounds[i] = start;
        }
    }
    return lowbounds;
}

int *linearSearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    for (int i = 0; i < m; i++) {
        lowbounds[i] = inline_linearSearch(0, n - 1, x, y[i]);
    }

    return lowbounds;
}

int *binarySearch(double *restrict x, int n, double *restrict y, int m) {
    int *lowbounds = (int *)malloc(m * sizeof(int));

    for (int i = 0; i < m; i++) {
        lowbounds[i] = inline_binarySearch(0, n - 1, x, y[i]);
    }

    return lowbounds;
}
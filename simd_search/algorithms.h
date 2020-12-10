#pragma once

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * creates a hash table based on the log_B() for each element in x
 * the base B is calculated such that the full range of a chosen
 * hashtable size is utilized
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *logHashSearch(double *restrict x, int n, double *restrict y, int m);

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * creates a list of indices such that the distance between the indicies
 * is equal in bytes to the size of a cachline, then uses this list to
 * search for targets
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *skiplistSearch(double *restrict x, int n, double *restrict y, int m);

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * creates a hash table based on the log_2() for each element in x
 * the base-2 exponent is extract directly from the exponent bits
 * of the values hashed
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *expHashSearch(double *restrict x, int n, double *restrict y, int m);

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * incrementally adjusts the bounds of each search so that correlations
 * between targets are utilized
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *huntLocateSearch(double *restrict x, int n, double *restrict y, int m);

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *linearSearch(double *restrict x, int n, double *restrict y, int m);

/**
 * @brief Performs a lower-bounds search for values from y in x.
 * 
 * @param x The data array; array to search for lower bounds within
 * @param n Number of elements in x
 * @param y The targets array; array of elements to find lower bounds for in x
 * @param m Number of elements in y
 * @return int* Array of lower bound indices, size is m
 */
int *binarySearch(double *restrict x, int n, double *restrict y, int m);
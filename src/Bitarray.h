#pragma once

#include <Rcpp.h>
#include <climits>
using namespace Rcpp;

#define BITS_IN_INT (sizeof(int) * 8)

/* Hat-tip to SO for these macros:
* "How to define and work with an array of bits in C?"
* https://stackoverflow.com/a/30590727/919872
*/

#define SetBit(A,k)     ( A[(k / BITS_IN_INT)] |= (1 << (k % BITS_IN_INT)) )
#define ClearBit(A,k)   ( A[(k / BITS_IN_INT)] &= ~(1 << (k % BITS_IN_INT)) )
#define TestBit(A,k)    ( A[(k / BITS_IN_INT)] & (1 << (k % BITS_IN_INT)) )

/*
 * Bitarray structure for storing binary data as array of ints
 * \param int** data an array of array of ints
 * \param nrow int storing number of elements to store -- usually a "row" of data from a sparse matrix
 * \param ncol int storing number of bytes needed to store all of the bits
 * \param nbits int storing number of individual bits to represent
 */
typedef struct  {
  int** data;
  int nrow, ncol, nbits;
} Bitarray;



int at(Bitarray x, int row, int bit);

void set(Bitarray x, int row, int bit);

int* get_row(Bitarray x, int row);

int* get_col(Bitarray x, int col);

/*
 * Helper method returning number of bytes need to store nbits for each machine
 * Calculated as the # of bits needed / # bits in an int
 */
int bytes_needed(int nbits);

/*
 * Allocate and zero out the space needed to store requested elements
 */
int** allocate_bytes(int nrow, int nbits);

/*
 * Constructor initializing a Bitarray struct and allocating spaces
 */
Bitarray new_bitarray(int nrow, int nbits);


/*
 * Destructor freeing memory allocated for Bitarray
 */
void free_bitarray(Bitarray x);


/* 
 * Find the mode of a bitarray
 * 
 * /param x A Bitarray for which to find the mode
 * /return Returns a bitarray with one row where each bit is set if more
 *         than 50% of the corresponding bits in x are set.
*/
Bitarray bitarray_mode(Bitarray x);

/*
 * convert bitarray to ngCMatrix to return to R
 */
Rcpp::S4 bitarray_to_ngCMatrix(Bitarray x);


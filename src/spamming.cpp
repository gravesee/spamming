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

// Lengths MUST be the same
int hamming_dist(int* a, int* b, size_t nels) {

  int dist = 0;
  for (int i = 0; i < nels; i++) {
    dist += __builtin_popcount(a[i] ^ b[i]);
  }

  return dist;
}

int bytes_needed(int size) {
  return size / BITS_IN_INT + 1;
}

// distance between
typedef struct  {
  int** data;
  int nrow, ncol;
} Bitarray;


void free_bitarray(Bitarray x) {
  for (int i = 0; i < x.nrow; i++) {
    free(x.data[i]);
  }
  free(x.data);
}


/* Translate a pattern matrix to an array of bits stored as ints
 * @param obj - a matrix of type ngCMatrix or lgCMatrix
 * @return a Bitarray struct
 */
Bitarray ngCMatrix_to_array(Rcpp::S4 obj) {

  // transposed matrix passed in so get the length of a row (which is a column here)
  int* dim = INTEGER(obj.slot("Dim"));

  int* i = INTEGER(obj.slot("i"));
  int* p = INTEGER(obj.slot("p"));

  // allocate all of the memory needed to store the bitsets
  int nrow = dim[1];
  int ncol = bytes_needed(dim[0]);
  int** m = (int**) calloc(nrow, sizeof(int*));

  // allocate arrays of ints representing each row
  for (int r = 0; r < nrow; r++) {
    m[r] = (int*) calloc(ncol, sizeof(int));
  }

  // loop over pointer indices
  for (int r = 0; r < nrow; r++) {

    int first_pos = p[r];
    int last_pos = p[r + 1];

    // Set the bits
    for (int c = first_pos; c < last_pos; c++) {
      SetBit(m[r], i[c]);
    }
  }

  Bitarray x;
  x.data = m;
  x.nrow = nrow;
  x.ncol = ncol;

  return x;

};


// [[Rcpp::export]]
IntegerMatrix hamming_ngCMatrix_x_only(Rcpp::S4 obj) {

  Bitarray x = ngCMatrix_to_array(obj);

  IntegerMatrix dist(x.nrow, x.nrow);

  for (int r = 0; r < x.nrow - 1; r++) {
    for (int nr = r + 1; nr < x.nrow; nr++) {
      dist(r, nr) = hamming_dist(x.data[r], x.data[nr], x.ncol);
    }
  }

  free_bitarray(x);

  return dist;

}


// [[Rcpp::export]]
IntegerMatrix hamming_ngCMatrix_x_and_y(Rcpp::S4 objx, Rcpp::S4 objy) {

  Bitarray x = ngCMatrix_to_array(objx);
  Bitarray y = ngCMatrix_to_array(objy);

  IntegerMatrix dist(x.nrow, y.nrow);

  for (int xr = 0; xr < x.nrow; xr++) {
    for (int yr = 0; yr < y.nrow; yr++) {

      dist(xr, yr) = hamming_dist(x.data[xr], y.data[yr], x.ncol);

    }
  }

  free_bitarray(x);
  free_bitarray(y);

  return dist;
}

#include "Bitarray.h"

Bitarray new_bitarray(int nrow, int nbits) {
  Bitarray x;
  x.data = allocate_bytes(nrow, nbits);
  x.nrow = nrow;
  x.nbits = nbits;
  x.ncol = bytes_needed(nbits);
  return x;
}

int bytes_needed(int nbits) {
  return nbits / BITS_IN_INT + 1;
}

void free_bitarray(Bitarray x) {
  for (int i = 0; i < x.nrow; i++) {
    free(x.data[i]);
  }
  free(x.data);
}

int** allocate_bytes(int nrow, int nbits) {
  int ncol = bytes_needed(nbits);
  int** m = (int**) calloc(nrow, sizeof(int*));
  
  // allocate arrays of ints representing each row
  for (int r = 0; r < nrow; r++) {
    m[r] = (int*) calloc(ncol, sizeof(int));
  }
  
  return m;
}

/* find the mode of a bitarray
*/
Bitarray bitarray_mode(Bitarray x) {
  
  Bitarray m = new_bitarray(1, x.nbits);
  
  // loop over rows and count bits
  for (int bit = 0; bit < x.nbits; bit++) {
    float cnt = 0.0;
    for (int row = 0; row < x.nrow; row++) {
      cnt += TestBit(x.data[row], bit) != 0;
    }
    // if cnt exceeeds 50% set the bit
    if ((cnt / x.nrow) > 0.50) {
      SetBit(m.data[0], bit);
    }
  }
  
  return m;
}

Rcpp::S4 bitarray_to_ngCMatrix(Bitarray x) {
  
  IntegerVector dim = IntegerVector::create(x.nrow, x.nbits);
  IntegerVector i;
  IntegerVector p;
  
  // loop over each bit
  int colptr = 0;
  for (int bit = 0; bit < x.nbits; bit++) {
    p.push_back(colptr);
    
    // loop over each row
    for (int r = 0; r  < x.nrow; r++) {
      // test if bit is set
      if (TestBit(x.data[r], bit) != 0) {
        i.push_back(r);
        colptr++;
      }
    }
  }
  p.push_back(colptr);
  
  // Create S4 Sparse Matrix
  S4 m("ngCMatrix");
  m.slot("i")  = i;
  m.slot("p") = p;
  m.slot("Dim") = dim;
  m.slot("Dimnames") = Rcpp::List(2);
  m.slot("factors") = Rcpp::List();
  m.attr("package") = "Matrix";
  
  return(m);
};



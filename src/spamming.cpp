#include <Rcpp.h>
#include <climits>
#include "Bitarray.h"
#include "BMM.h"
//#include "ExpectationMaximization.h"
using namespace Rcpp;

// Lengths MUST be the same
int hamming_dist(int* a, int* b, size_t nels) {

  int dist = 0;
  for (int i = 0; i < nels; i++) {
    dist += __builtin_popcount(a[i] ^ b[i]);
  }

  return dist;
}

/* Translate a pattern matrix to an array of bits stored as ints
 * @param obj - a matrix of type ngCMatrix or lgCMatrix
 * @return a Bitarray struct
 */
// [[Rcpp::export]]
SEXP ngCMatrix_to_array_test(Rcpp::S4 obj) {

  // transposed matrix passed in so get the length of a row (which is a column here)
  
  int* dim = INTEGER(obj.slot("Dim"));
  int* i = INTEGER(obj.slot("i"));
  int* p = INTEGER(obj.slot("p"));

  // allocate all of the memory needed to store the bitsets
  //Bitarray x = new_bitarray(dim[1], dim[0]);
  Bitarray x = new_bitarray(dim[0], dim[1]);
  
  // loop over pointer indices
  for (int bit = 0; bit < x.nbits; bit++) {
    int start = p[bit];
    int stop = p[bit + 1];
    // loop over rows
    for (int r = start; r < stop; r++) {
      int row = i[r];
      set(x, row, bit);
    }
  }
  
  free_bitarray(x);
  return R_NilValue;
};


Bitarray ngCMatrix_to_array(Rcpp::S4 obj) {

  int* dim = INTEGER(obj.slot("Dim"));
  int* i = INTEGER(obj.slot("i"));
  int* p = INTEGER(obj.slot("p"));
  
  Bitarray x = new_bitarray(dim[0], dim[1]);
  
  // loop over pointer indices
  for (int bit = 0; bit < x.nbits; bit++) {
    int start = p[bit];
    int stop = p[bit + 1];
    // loop over rows
    for (int r = start; r < stop; r++) {
      int row = i[r];
      set(x, row, bit);
      //SetBit(x.data[row], bit);
    }
  }
  
  return x;
  
};

// [[Rcpp::export]]
Rcpp::IntegerVector hamming_ngCMatrix_x_only(Rcpp::S4 obj) {
  
  Bitarray x = ngCMatrix_to_array(obj);
  
  IntegerVector dist((x.nrow * (x.nrow - 1)) / 2);
  
  int i = 0;
  for (int r = 0; r < x.nrow - 1; r++) {
    for (int nr = r + 1; nr < x.nrow; nr++, i++) {
      dist[i] = hamming_dist(x.data[r], x.data[nr], x.ncol);
    }
  }
  
  free_bitarray(x);
  
  // Create an S3 dist object
  dist.attr("Size")  = IntegerVector::create(x.nrow);
  dist.attr("call")  = R_NilValue;
  dist.attr("class") = CharacterVector::create("dist");
  dist.attr("Diag")  = LogicalVector::create(false);
  dist.attr("Upper") = LogicalVector::create(false);
  
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



// [[Rcpp::export]]
S4 hamming_find_mode(Rcpp::S4 obj) {
  
  Bitarray x = ngCMatrix_to_array(obj);
  Bitarray m = bitarray_mode(x);
  
  S4 out = bitarray_to_ngCMatrix(m);
  
  free_bitarray(x);
  free_bitarray(m);
  
  return out;
}

// [[Rcpp::export]]
S4 test_conversion(Rcpp::S4 obj) {
  Bitarray x = ngCMatrix_to_array(obj);
  S4 out = bitarray_to_ngCMatrix(x);
  free_bitarray(x);
  return(out);
}

// IntegerVector BMM(Rcpp::S4 obj, int K, int repetitions) {
//   Bitarray x = ngCMatrix_to_array(obj);
//   
//   Rprintf("Creating EM object\n");
//   ExpectationMaximization em = new_ExpectationMaximization(x, K);
//   
//   Rprintf("EM Step\n");
//   PerformEMstep(&em, repetitions);
//   
//   Rprintf("Getting cluster assignments\n");
//   IntegerVector result(x.nrow);
//   for (int r = 0; r < x.nrow; r++) {
//     result(r) = GetCluster(&em, x, r);
//   }
//   
//   Rprintf("Free EM\n");
//   free_ExpectationMaximization(em);
//   
//   Rprintf("Free Bitarray\n");
//   free_bitarray(x);
//   
//   return result;
// }





// [[Rcpp::export]]
SEXP BMM_v2(Rcpp::S4 obj, int K, int max_iter, int verbose) {
  Bitarray x = ngCMatrix_to_array(obj);
  
  BMM_Result res = em(x, K, max_iter, verbose);
  
  // load stuff into R object
  Rcpp::List protos;
  Rcpp::NumericVector pis(K);
  for (int k = 0; k < K; k++) {
    Rcpp::NumericVector proto(x.nbits);
    
    for (int i = 0; i < x.nbits; i++) {
      proto[i] = res.protos[k][i];
    }
    
    protos.push_back(proto);
    
    pis[k] = res.pis[k];
  }
  
  Rcpp::IntegerVector cluster(x.nrow);
  for (int n = 0; n < x.nrow; n++) {
    cluster[n] = res.cluster[n];
  } 

  free_BMM_Result(&res);
  
  return Rcpp::List::create(
    Rcpp::Named("prototypes") = protos,
    Rcpp::Named("pis") = pis,
    Rcpp::Named("cluster") = cluster,
    Rcpp::Named("ll") = res.ll
  );
  
};




#include "ExpectationMaximization.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <limits.h>
// using namespace Rcpp;

ExpectationMaximization new_ExpectationMaximization (Bitarray x, int K) {
  
  ExpectationMaximization em;
  
  // Set the Bitarray as the data to cluster
  em.x = x;
  
  // K (number of clusters)
  // Initialize D (dimensions)
  // N (number of items to be clustered)
  em.K = K;
  em.D = x.nbits;
  em.N = x.nrow;
  
  // Initialize pi_1, ..., pi_K (mixing coefficients) to uniform values
  em._pi = (double*) malloc(em.K * sizeof(double));
  for (int k = 0; k < em.K; k++) {
    em._pi[k] = 1.0/em.K;
  }
  
  // Initialize mu_1, ..., mu_K (cluster's bit distributions)
  // allocate memory
  em._mu = (double**) calloc(em.K, sizeof(double*));
  
  for (int k = 0; k < em.K; k++) {
    em._mu[k] = (double*) calloc(em.D, sizeof(double));
    
    double normalizationFactor = 0.0;
    for (int i = 0; i < em.D; i++) {
      
      GetRNGstate();
      em._mu[k][i] = (unif_rand() * 0.5) + 0.25;
      PutRNGstate();
      
      normalizationFactor += em._mu[k][i];
    }
    
    for (int i = 0; i < em.D; i++) {
      em._mu[k][i] /= normalizationFactor;
    }
  }
  
  // Initialize z_1, ..., z_N (latent variable over cluster membership of x_1, ... x_N)
  em._z = (double**) calloc(em.N, sizeof(double*));
  for (int n = 0; n < em.N; n++) {
    em._z[n] = (double*) calloc(em.K, sizeof(double));
  }
  
  return em;
};

void free_ExpectationMaximization(ExpectationMaximization em) {
  
  // Need to free _pi, _mu and _z;
  // _pi
  free(em._pi);
  
  // _mu
  for (int k = 0; k < em.K; k++) {
    free(em._mu[k]);
  }
  free(em._mu);
  
  // _z
  for (int n = 0; n < em.N; n++) {
    free(em._z[n]);
  }
  free(em._z);
  
};

void PerformEMstep(ExpectationMaximization* em, int repetitions) {
  
    for (int i = 0; i < repetitions; i++) {
      ExpectationStep(em);
      MaximizationStep(em);
    }
  
};

void ExpectationStep(ExpectationMaximization* em) {
  
    for (int n = 0; n < em->N; n++)
    {
      double normalizationFactor = 0.0;
      
      for (int k = 0; k < em->K; k++) {
        
        em->_z[n][k] = ExpectationSubstep(em, n, k);
        normalizationFactor += em->_z[n][k];
      }
      
      // Re-normalize z[n, k]
      for (int k = 0; k < em->K; k++)
      {
        if (normalizationFactor > 0.0) {
          em->_z[n][k] /= normalizationFactor;
        }
        else {
          em->_z[n][k] = 1.0 / (float) em->K;
        }
      }
    }
    
  
};

void MaximizationStep(ExpectationMaximization* em) {
  
  // Update pi_1, ..., pi_k
  for (int k = 0; k < em->K; k++)
  {
    em->_pi[k] = Nm(em, k) / (double) em->N;
  }
  
  // Update mu_1, ..., mu_K
  for (int k = 0; k < em->K; k++) {
    double* averageX_k = AverageX(em, k);
    
    for (int i = 0; i < em->D; i++)
    {
      em->_mu[k][i] = averageX_k[i];
    }
  }
  
};

double ExpectationSubstep(ExpectationMaximization* em, int n, int k) {
  double z_nk = em->_pi[k];
  for (int i = 0; i < em->D; i++) {

    int bit_is_set = (TestBit(em->x.data[n], i)!=0);
    //z_nk *= pow(em->_mu[k][i], bit_is_set) * pow(1.0 - em->_mu[k][i], 1.0 - bit_is_set);
    
    z_nk += log(pow(em->_mu[k][i], bit_is_set)) + log(pow(1.0 - em->_mu[k][i], 1.0 - bit_is_set));
  
  }
  
  return exp(z_nk);
};

double* AverageX(ExpectationMaximization* em, int m) {
  
  double* result = (double*) calloc(em->D, sizeof(double));
  
  for (int i = 0; i < em->D; i++){
    for (int n = 0; n < em->N; n++) {
      
      result[i] += em->_z[n][m] * (TestBit(em->x.data[n], i)!=0);
    }
  }
  
  double currentNm = Nm(em, m);
  for (int i = 0; i < em->D; i++) {
    result[i] /= currentNm;
  }
  
  return result;
  
};

double Nm(ExpectationMaximization* em, int m) {
  
  double result = 0.0;
  
  for (int n = 0; n < em->N; n++) {
    result += em->_z[n][m];
  }
  
  return result;
  
};

// em object, bitarray object, row from bitarray to get cluster id 
int GetCluster(ExpectationMaximization* em, Bitarray x, int n) {
  
  double maxClusterSum = DBL_MIN;
  int maxCluster = -1;
  
  for (int k = 0; k < em->K; k++) {
    
    double currentClusterSum = 0.0;
    for (int i = 0; i < em->D; i++) {
      currentClusterSum += (TestBit(x.data[n], i)) ? em->_mu[k][i] : 1.0 - em->_mu[k][i];
    }
    
    if (currentClusterSum > maxClusterSum) {
      maxClusterSum = currentClusterSum;
      maxCluster = k;
    }
  }
  
  return maxCluster;
};
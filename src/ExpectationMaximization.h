#include "Bitarray.h"

// porting C# code found here: https://github.com/manfredzab/bernoulli-mixture-models

typedef struct {
  
  Bitarray x;
  int K, D, N;
  double* _pi;
  double** _mu;
  double** _z;
  
} ExpectationMaximization; // Expectation Maximation struct

ExpectationMaximization new_ExpectationMaximization (Bitarray x, int K);

void free_ExpectationMaximization(ExpectationMaximization em);

void PerformEMstep(ExpectationMaximization* em, int repetitions);

void ExpectationStep(ExpectationMaximization* em);

void MaximizationStep(ExpectationMaximization* em);

double ExpectationSubstep(ExpectationMaximization* em, int n, int k);

double* AverageX(ExpectationMaximization* em, int m);

double Nm(ExpectationMaximization* em, int m);

int GetCluster(ExpectationMaximization* em, Bitarray x, int n);
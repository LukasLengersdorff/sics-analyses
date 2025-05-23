data {
  int<lower=0> n;        // number of dimensions
  vector[n] mu;          //
  matrix[n,n] Sigma;    //
  //int sample_posterior;  //

  /*vvv SICS vvv*/
  // A: equality constraints
  int<lower=0> SICS_mA;
  matrix[n,SICS_mA] SICS_MA;
  vector[SICS_mA] SICS_a;
  // B: inequality constraints
  int<lower=0> SICS_mB;
  matrix[n,SICS_mB] SICS_MB;
  vector[SICS_mB] SICS_b;
  // C: interval constraints
  int<lower=0> SICS_mC;
  matrix[n,SICS_mC] SICS_MC;
  vector[SICS_mC] SICS_c0;
  vector[SICS_mC] SICS_c1;
  // U: unconstrained part
  int<lower=0> SICS_mU;
  matrix[n,SICS_mU] SICS_MU;
  /*^^^ SICS ^^^*/
}

transformed data{
  /*vvv SICS vvv*/
  vector[n]  SICS_MA_x_a;
  if (SICS_mA > 0) {
    SICS_MA_x_a = SICS_MA * SICS_a;
  }
  /*^^^ SICS ^^^*/
}

parameters {
  /*vvv SICS vvv*/
  vector<lower=SICS_b>[SICS_mB] SICS_wB;
  vector<lower=SICS_c0, upper=SICS_c1>[SICS_mC] SICS_wC;
  vector[SICS_mU] SICS_wU;
  /*^^^ SICS ^^^*/
}

transformed parameters {
  vector[n] theta;

  /*vvv SICS vvv*/
  {vector[n] SICS_wsum = rep_vector(0, n);
   if (SICS_mA > 0) SICS_wsum += SICS_MA_x_a;
   if (SICS_mB > 0) SICS_wsum += SICS_MB * SICS_wB;
   if (SICS_mC > 0) SICS_wsum += SICS_MC * SICS_wC;
   if (SICS_mU > 0) SICS_wsum += SICS_MU * SICS_wU;
   theta = SICS_wsum;
  }
  /*^^^ SICS ^^^*/
}

model {
  theta ~ multi_normal(mu, Sigma);
}

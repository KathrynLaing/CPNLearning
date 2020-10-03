#include "10_DirichletSample.h"
using namespace std;


std::vector<double> DirSample(std::vector<double> beta, int k, int N){
  //This Function takes k beta values and draws N times from the associated Dir distribution
  //sample is output as a N*k vector listing each sample in order
  double InputBeta[k],sample[k];
  //We store the total sample of N draws (each length k) in totalSample
  std::vector<double> totalSample(N*k);
  //InputBeta is the vector of input parameters
  for(int i=0;i<k;i++){
    InputBeta[i]=beta[i];
  }
  //Set up the random generator, using gsl's dirichlet random generator
  gsl_rng_env_setup();
  gsl_rng *r;
  r=gsl_rng_alloc(gsl_rng_mt19937);
  struct timeval tv;
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  gsl_rng_set(r, mySeed);
  //Repeat N times:
  for(int i=0;i<N;i++){
    //Randomly draw from the Dirichlet distribution of order k with parameters InputBeta.
    //Store the draw in sample
    gsl_ran_dirichlet(r,k,InputBeta,sample);
    //Add this draw to the end of totalSample
    for(int j=0;j<k;j++){
      totalSample[(k*i +j)]=sample[j];
    }
  }
  gsl_rng_free(r);
  //totalSample now lists the N samples in order
  return totalSample;
}

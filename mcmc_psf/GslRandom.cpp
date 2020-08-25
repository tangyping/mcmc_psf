#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
using namespace std;

#include "GslRandom.h"

//constructor
GslRandom::GslRandom(unsigned seed)
{
  //create the generator
  r = gsl_rng_alloc(gsl_rng_mt19937);  //using the MT19937 generator

  //gsl_rng * g = gsl_rng_alloc (gsl_rng_taus);

  //set the random number seed based on the current system time
  //which is seconds since a long time ago

  //set the generator
  gsl_rng_set(r,seed);
}

//return a gaussian deviate with variance=1
double GslRandom::gaussDeviate()
{
  return gsl_ran_ugaussian(r);
}

double GslRandom::betaDeviate(double a)
{
  return gsl_ran_beta(r, a, a);
}

//returns a uniform deviate on the range [a,b]
double GslRandom::uniformDeviate(double a, double b)
{
  return gsl_ran_flat(r, a, b);
}

double GslRandom::expDeviate(double mu)
{
  return gsl_ran_exponential(r, mu);
}

GslRandom::~GslRandom()
{
  gsl_rng_free(r);
}

  

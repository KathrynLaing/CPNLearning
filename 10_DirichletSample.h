#ifndef DIRICHLETSAMPLE_H
#define DIRICHLETSAMPLE_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
using namespace std;


std::vector<double> DirSample(std::vector<double> alpha, int k, int N);

#endif

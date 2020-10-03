#ifndef CPNLTEST_H
#define CPNLTEST_H

#include "10_DirichletSample.h"
#include "cpnInput.h"
#include "cpnOutput.h"
#include "cpnlearn.h"
#include "domtest.h"
#include <random>
#include <chrono>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

double dataAgreementWithFlips(string filename, int D, cpn cpnet);

double dataOrderCompatibleWithCPN(string filename, int D, cpn cpnet, long long int nTests);

#endif
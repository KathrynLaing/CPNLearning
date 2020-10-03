#ifndef VALIDATIONMEASURES_H
#define VALIDATIONMEASURES_H

#include "10_DirichletSample.h"
#include "cpnInput.h"
#include "cpnOutput.h"
#include "cpnlearn.h"
#include "domtest.h"
#include <random>
#include <chrono>
#include <sys/time.h>
using namespace std;

vector<double> EntailmentAgreement(cpn NLearned, vector<vector<int>> DQ1, vector<vector<int>> DQ2);

int PrefGraphMisorientation(cpn NOriginal, cpn NLearned);

#endif
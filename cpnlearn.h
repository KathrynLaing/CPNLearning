#ifndef CPNLEARN_H
#define CPNLEARN_H

#include <math.h>
using namespace std;

LearningOutput cpnLearn(int nvar, std::vector<int> data, std::vector<double> beta, int DataSize, int mcReps, double signif, bool StartEmpty=true, vector<int> StartStructure={0}, string CpnSeqFilename="/file/path/output.csv");

#endif

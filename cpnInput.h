#ifndef CPNINPUT_H
#define CPNINPUT_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

struct cpn
{
    int size;
    std::vector<int> structure;
    std::vector<int> domains;
    std::vector<int> cptEntries;
    std::vector<int> cptBreaks;
};

struct LearningOutput
{
    vector<cpn> LearnedCpnSeq;
    vector<double> FinalScores;
    double TimeElapsed;
};

struct cpn readCpn(string filename, int n);

vector<cpn> readMultCpn(string filename, int n);

#endif

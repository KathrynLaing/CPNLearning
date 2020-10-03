#ifndef DOMTEST_H
#define DOMTEST_H

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <limits>

using namespace std;

struct List
{
    bool result;
    long long int OT;
};

vector<int> CPTRow(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, int x, vector<int> PA);

double Rank(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, vector<int> o);

std::vector<double> RankDifferences(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK);

List RSDQRankPriority(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, vector<int> o1, vector<int> o2);

#endif

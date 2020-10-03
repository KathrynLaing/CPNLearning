# CP-Net Learning

This repository contains the code for my CP-net learning method (and associated functions), the details of which can be found in my thesis:

*Laing, K. (2020). Conditional Preference Networks: Efficient Dominance Testing and Learning. PhD Thesis. School of Mathematics, University of Leeds, UK.*

This CP-net learning method takes a set of observed user choices and learns a conditional preference network (CP-net) model representation of user preference. Note that, in its current form we can only use this method on problems with binary variables.

The C++ files contain the necessary functions to perform learning and to test the outputs of our learning. The numbering of these files indicates the hieirarchy of the code. That is, file 20 can be dependent upon functions from file 10 but not vice versa.

CompilingInstructions.txt gives the necessary commands to compile these scripts. Note that the gsl library is required.

Below we will outline the functions given by each script and any additional information required to utilise them. We will also discuss the formatting methods for objects such as outcomes, CP-nets, and data. Finally, we will discuss...

## Functions

**10_DirichletSample.cpp**

This file contains the function for randomly drawing samples from a given Dirichlet distribution.
```
DirSample(beta, k, N)
```
`Beta` is the vector of Dirichlet parameters. `k` is the number of parameters (|beta|), that is, the order of the Dirichlet distribution. `N` is the number of draws we want.

We use the `gsl_ran_dirichlet` function from the gsl library to randomly generate these samples.

Each sample is a vector of length k. We output the N samples as a single vector of length N*k with each sample listed end to end. That is, entries 0 to k-1 is draw 1 from the Dirichlet distribution. Entries k to 2k-1 constitute the second draw and so on.

**25_CPNImport.cpp**

This file contains two functions for reading in CP-nets from text files.
```
readCpn(filename, n)
readMultCpn(filename, n)
```
`filename` is the location of the text file. `n` is the number of variables in the CP-nets.

`readCpn` takes a text file and outputs a cp-net. `readMultCpn` takes a text file that may have more than one CP-net and outputs a vector of the CP-nets. 

The required format of text files is explained in the format section below

These functions, in combination with the following functions which write CP-nets to file, allow us to store and call CP-nets in a convenient form as needed.

**26_CPNOutput.cpp**

This file contains a function for writing CP-nets to file.
```
writeCpn(filename, N)
```
`filename` is the location of the text file we write to. `N` is the CP-net of interest. 

If the exists already, then the CP-net is appended to the end of the file. The ouput format is the same as the input format for `readCpn`, given at the bottom in detail. We use the same format so that all saved CP-nets are in the same form and can all be read in by the same function (`readCpn` or `readMultCpn`).

**30_cpnLearning.cpp** 

This file contains my CP-net learning function
```
LearningOutput cpnLearn(n, data, beta, DataSize, mcReps, alpha, StartEmpty=true, StartStructure={0}, CpnSeqFilename="/file/path/output.csv")
```
-`n` Number of variables
-`data` Observed user choices data (format given below)
-`beta` Vector of Dirichlet prior parameters
-`mcReps` Number of dirichlet samples used in score estimation
-`alpha` Change threshold
-`StartEmpty` Boolean variable, true if and only if we want to learn from the empty structure
-`StartStructure` The specified starting structure if StartEmpty=false (set to {0} otherwise)
-`CpnSeqFilename` The location we want the learning results to be written to

This function takes the observed data and input parameters and applies our learning procedure. The output structure contains a sequence of the CP-nets obtained through learning, a vector of n variable scores, and the learning time elapsed. 

The sequence of CP-nets is the CP-nets obtained after each edge change performed, followed by the final learned CP-net. In the case where at least one edge change is performed, this means the final CP-net is given twice.

The variable score vector gives the table scores of each variable in the final learned structure. The product of this vector gives the CP-net score of the learned structure

These output are all also written to CpnSeqFilename.

**35_DominanceTesting.cpp**
different lex encoding cause 1/2 not 0/1
precisions
 (Note that the R function for this method can be found in ... repository and the other DT fns can be found....).

**40_cpnlTest.cpp**

This file contains the functions for evaluating DFA and DOC scores between a CP-net and a given test set.
```
dataAgreementWithFlips(filename, D, cpnet)
```
This calculates the DFA between the CP-net `cpnet` and the data located at `filename`. `D` is the number of data points in our test set. The format of the data is described below.
```
dataOrderCompatibleWithCPNfilename, D, cpnet, nTests)
```
This function calculates the DOC between the CP-net `cpnet` and the data located at `filename`. `D` is the number of data points in our test set. The data is assumed to be formatted in the same way as for DFA. `nTests` is the number of outcome pairs we want to test to calculate DOC. This cannot be more than 2^(n-1)\*(2^n-1). If a larger value is entered we convert to the maximal value 2^(n-1)\*(2^n-1). 

This function requires dominance testing. We use rank pruning and suffix fixing. The function itself (in C++) is given in 35_DominanceTesting.cpp. Using this function requires a slight adjustment. All of our CP-nets in C++ have breaks vector starting from 0. Our previous R encoding of CP-nets had breaks vector starting at 1. Due to an issue in translation from R to C++ of our dominance testing function, we must add 1 to all entries of `breaks` before passing the CP-net to the dominance testing function.

**50_CPNLFurtherValidation** contains two testing functions- entailment agreement and preference graph similarity

EntailmentAgreement(N_L, DQ1, DQ2)

Takes the learned CP-net and a set of N preferences (o1>o2). DQ1 is the list of the preferred outcomes (o1) and DQ2 is the list of the not-preferred outcomes (o2) in the same order. These preferences are all entailed by N_T.

The output is three integers, specifying how many of these preferences are 
a)entailed by N_L
b)contradicted by N_L (the opposite preference is entailed)
c)neither entailed nor contradicted by N_L

[Note these 3 options are mutually exclusive]

Determining this requires dominance testing. We use rank pruning and suffix fixing. The function itself (in C++) is given in.... (Note that the R function for this method can be found in ... repository and the other DT fns can be found....). Using this function requires a slight adjustment. All of our CP-nets in C++ have breaks vector starting from 0. R CP-nets have breaks vector starting at 1. Due to an issue in translation from R to C++, we must add 1 to all entries of `breaks' before passing the CP-net to the dominance testing function.

In order to obtain the test set, we simply generate random pairs of outcomes (making usre to never generate the same pair twice) and determine whether o1>o2 or o2>o1 is entailed by N_T. If a direction is entailed, add this to the set of test preferences. Continue this process until your test preference set is as large as necessary. For details on the limitations of sizes of such sets, see my thesis.

In order to determine if a proposed preference is entailed by N_T, we use the same dominance testing function used in EntailmentAgreement, this will include the requires adjustment of breaks.

PrefGraphMisorientation(N_T, N_L)

Takes N_T and N_L as cpn structures and returns the number of preference graph edges that are oriented differently for N_T and N_L (divide this by 2^(n-1)(2^n-1) to get PG similarity. 

This function uses the the lexicographic encoding of parent assignments heavily. Particularly in finding the right entry of CPTs.

To understand why this algorithm works, we need to show the following result. The number of X-flips that are oriented differently in the preference graph of N_T and N_L is as follows. Let P_T and P_L be the respective parent sets of X. Let I denote their intersections and O the variables not in either parent set.

Sum_I 2^|O|* ((#1>2 pref rules in N_T * #2>1 pref rules in N_L) + (#2>1 pref rules in N_T * #1>2 pref rules in N_L))

Sum_I - summing over all possible intersection value assignments
#1>2 pref rules in N_T - the number of rows in CPT(X) in N_T such that I=i (the specific assignment we are considering) and the preference rule is 1>2

This is true because O is not in either parental assignment and so does not affect orientation. In both preference graph edge orientation is the same for all 2^|O| assignments of O so we multiply our answer by 2^|O| to account for all possible O assignments. Given now that I=i, the orientation of an X-flip in N_T depends only on P_T\I. Given P_T\I =p the X flips will be in the same direction for all assignments to P_L\I. If the CPT(X) rule for P_T\I =p & I=i is 1>2 then these edge orientations will not match N_G for exactly the P_L\I assignments, q, such that P_L\I=q & I=i implies the rule 2>1 in N_L. Thus, to find the number of instances of misoriented edges we have to multiply the number of 1>2 cases of P_T\I by the number of 2>1 cases of P_L\I and also the number of 2>1 cases of P_T\I by the number of 1>2 cases of P_L\I.

This is the fundamental result used in this function. However in some cases it is simpler. e.g. when there is no intersection or P_L=P_T. Thus, we separate into cases.


Format
Outcomes are each an integer between 0 and 2^n-1. Data is a vector of integers representing the outcomes we have observed being chosen.

The number of an outcome is the lexicographic rank of it viewed as a vector. Binay variables mean that, once a variable order is fixed, every outcome is a 0/1 vector. Then outcome (0,0,0,0) is 0, (0,0,0,1) is 1, (0,0,1,0) is 2, and so on. We assume that the beta values are given in this lexicographic order and so are the draws from the Dirichlet distribution. Note that the value of an outcome (given its vector form/ the variable assignments) can be calculated via SUM o[i]2^(n-i). We use this correspondence a lot throughout to determine which entries of samples we should be considering.

All structures are given as n*n length vectors. This is essentially the adjacency matrix where rows are given end-to-end instead of on top of one another

A CPT is a 0/1 vector of length 2^|Pa|. Consider the m parents. Each can take value 0 or 1, so parent assignments are lenght m 0/1 vectors (their ordering is enforced by the fixed variable ordering). These vectors are ordered lexicographically. The ith entry of CPT vector is the preference rule under the parent assignment that is ith in the lexicographic order. A 0 entry in the CPT implies the preference rule 1>2 and an entry of 1 implies 2>1 preference (Note that every variable is considered to have domain {1,2} in our code as value labels are irrelevant). This representation is used in our code as well as the entries and breaks representation (which encodes ALL CPTs for the whole CP-net).


## Format
In both cases all CP-nets must be written in the files as an adjacency matrix (in nxn form), then the domain sizes vector on the next line, then the CPT entries vector on the next line, then the CPT breaks vector on the next line. All in csv format. This format of a CP-net matches the one we used in .....rcpp repository.

cpnInput.h

This file defines two new data structures in C++. The first is cpn (CP-nets). We are here using the same format for CP-nets as in the ....rcpp repository.... A cpn consists of size (number of variables), structure (the adjacency matrix of the CP-net, now flattened to a vector with rows written one after another), domains (vector of variable domain sizes), cpt Entries, and cptBreaks (the latter two are defined as in ...rcpp repository.. only now the breaks vector starts at 0 instead of 0 to be in line with C++ indexing. That is, take one away from each entry of Rcpp break vectors.

This structure gives us a compact and convenient way to deal with CP-nets in C++ as they have so many separate parts. Throughout, when we talk about using CP-nets or reading/writing CP-nets from C++ we will be using this structure type.

This file is identical to our text files storing CP-nets in the R case. It gives Adjacency matrix, Domain vector, Entries vector, and then Breaks vector, each on new lines. Note that the breaks vector starts at 1, as in the previous case, to match up with R indexing from 1 convenience. However, in C++ the breaks vectors we use start at 0.

The second struct is LearningOutput. This is made in particular as an object to return from our learning algorithm. It has a vector of cpns, which means we can give the CP-nets visited by learning after each edge change (and then the final learned structure). It also has a vector of doubles to give the variable table scores for the learned structure (the product of this vector is the learned CP-net score). Finally, it has a double to give the time elapsed while learning. This struct gives us a single object to return from learning which still allows for convenient access to all of the relevant information about learning performance.

We expect this data to consist of a single comma separated line of D observed outcome choices. Each outcome choice is a single integer that indicates which outcome was chosen. Again we are assuming a lexicographic enumeration of outcomes with a fixed variable order that matches the variable order in the CP-net.

 Again we use a lexicographic enumeration of outcomes (and also outcome pairs) in the function, often translating between outcomes and their enumerations.
 
CPN structures - encoding in C++
cpt entries and breaks
CPT indexit, always 12 or 01?

cyclic update - > thesis

#include "40_cpnlTest.h"


double dataAgreementWithFlips(string filename, int D, cpn cpnet){
    //Calculates DFA between data and cpnet
    //filename tells us the location of the data - Assumed to be written on 1 comma separated line, each entry is an observed outcome
    //D tells us the number of data points
    //cpnet is the cpn of interest
    //first we read in the data and enter it into vector data
    vector<int> data(D);
    int line=1;
    int cols=D;
    int d;
    ifstream file(filename);
    char dummy;
    for (int i = 0; i < line; i++){
        for (int j = 0; j < cols; j++){
            file >> d;
            data[j]=d;
            if (j < (cols - 1)){
                file >> dummy;
            }
        }
    }
    file.close();
    //nvar is the number of variables in the CP-net
    int nvar=cpnet.size;
    //convert the data into counts - instead of a list of observed outcomes we have a vector where the ith entry is the number of times outcome i was observed
    //outcomes are assumed to be in lexicographic order
    vector<int> counts(pow(2,nvar),0);
    for(int i=0;i<D;i++){
        counts[data[i]]+=1;
    }
    //LexMult is the lexicographic multipliers of the variables (the mutlipliers we use to convert between outcomes and their enumeration)
    vector<int> LexMult(nvar);
    for(int i=0;i<nvar;i++){
        LexMult[i]=pow(2,nvar-i-1);
    }
    //AgreementScore will be the DFA score numerator
    int AgreementScore=0;
    //adjMat is the adjacency matrix of the CP-net structure. entries and breaks are the cpt details
    vector<int> adjMat=cpnet.structure;
    vector<int> entries=cpnet.cptEntries;
    vector<int> breaks=cpnet.cptBreaks;
    //for each variable, X, calculate (and add) the the data differences over X flips
    for(int i=0;i<nvar;i++){
        //parents - 0/1 vector giving the parents of X
        //nPa - number of parents
        vector<int> parents(nvar);
        int nPa=0;
        for(int j=0;j<nvar;j++){
            parents[j]=adjMat[j*nvar+i];
            nPa+=parents[j];
        }
        //Pa - an nPa vector giving the indices of the parents
        vector<int> Pa(nPa);
        int counter=0;
        for(int j=0;j<nvar;j++){
            if(parents[j]==1){
                Pa[counter]=j;
                counter+=1;
            }
        }
        // nother is the number of variables which are not X or parents of X
        int nOther=nvar-1-nPa;
        // otherWeights - list of all possible weights these other variables can add to the lexicographic posn of an outcome
        //otherAssts - cycles through all possible assignments to these other variables
        std::vector<int> otherWeights(pow(2,nOther));
        std::vector<int> otherAssts(nOther,0);
        std::vector<int> otherMult(nOther);
        if(nOther>0){
            int counter=0;
            for(int k=0;k<nvar;k++){
                if(k!=i){
                    if(parents[k]==0){
                        otherMult[counter]=LexMult[k];
                        counter+=1;
                    }
                }
            }
            for(int k=0;k<pow(2,nOther);k++){
                int weight=0;
                for(int l=0;l<nOther;l++){
                    weight+= otherAssts[l]*otherMult[l];
                }
                otherWeights[k]=weight;
                if(k!=(pow(2,nOther)-1)){
                    int max=0;
                    for(int l=0;l<nOther;l++){
                        if(otherAssts[l]==0){
                            if(l>max){
                                max=l;
                            }
                        }
                    }
                    otherAssts[max]=1;
                    for(int l=max+1;l<nOther;l++){
                        otherAssts[l]=0;
                    }
                }
            }
        }
        //For each parent assignment to Pa(X) we find the data differences over all X flips under this assignment
        //PaAsst cycles through all possible parental assignments
        vector<int> PaAsst(nPa,0);
        for(int j=0; j<pow(2,nPa);j++){
            //fixedweight is the weight contributed by parent assignment PaAsst to the lexicographic position of an outcome
            int fixedweight=0;
            for(int k=0;k<nPa;k++){
                fixedweight+=PaAsst[k]*LexMult[Pa[k]];
            }
            //group1/2 are the lists of the lexicographic positions of outcomes with Pa(x)=PaAsst and X=1/2
            vector<int> group1(pow(2,nOther));
            vector<int> group2(pow(2,nOther));
            if(nOther==0){
                group1[0]=fixedweight;
                group2[0]=fixedweight+LexMult[i];
            }
            else{
                for(int k=0;k<pow(2,nOther);k++){
                    group1[k]=fixedweight+otherWeights[k];
                    group2[k]=fixedweight+otherWeights[k]+LexMult[i];
                }
            }
            //value1/2 is the number of data observations of outcomes in group 1/2
            int value1=0;
            int value2=0;
            for(int k=0;k<pow(2,nOther);k++){
                value1+=counts[group1[k]];
                value2+=counts[group2[k]];
            }
            //if the CPT rule corresponding to PaAsst is 1>2, then the data difference is value1-value2
            // Add this to our score
            if(entries[breaks[i]+2*j]==1){
                AgreementScore+=value1-value2;
            }
            else{//if the CPT rule is 2>1, then the data difference is value2-value1
                AgreementScore+=value2-value1;
            }
            //move to next parent assignmnent
            if(j!=(pow(2,nPa)-1)){
                int max=0;
                for(int k=0;k<nPa;k++){
                    if(PaAsst[k]==0){
                        max=k;
                    }
                }
                PaAsst[max]=1;
                for(int k=max+1;k<nPa;k++){
                    PaAsst[k]=0;
                }
            }
            //We have added the data differences for all PaAsst X-flips to the score
        }
        //After cycling through all parental assignments, we have now added all X-flip data differences
    }
    //AgreementScore is now the sum of all variable flip data differences
    //To obtain DFA, we must divide this by n*#data points
    double ScaledFlipAgreement = (double) AgreementScore/ ((double) nvar* (double) D);
    return ScaledFlipAgreement;
}


double dataOrderCompatibleWithCPN(string filename, int D, cpn cpnet, long long int nTests){
    //Calculates DOC between the cpnet and data located at filename. nTests the level of approximation we will use in our calculation.
    //filename - location of data
    //D - number of data points in file
    //cpn - CP-net of interest
    //nTests - number of comparisons we are going to do (max possible number of tests is number of unordered pairs of outcomes)
    //first we read in the data and enter it into vector data
    vector<int> data(D);
    int line=1;
    int cols=D;
    int d;
    ifstream file(filename);
    char dummy;
    for (int i = 0; i < line; i++){
        for (int j = 0; j < cols; j++){
            file >> d;
            data[j]=d;
            if (j < (cols - 1)){
                file >> dummy;
            }
        }
    }
    file.close();
    //nvar is the number of variables in the CP-net
    int nvar=cpnet.size;
    //convert the data into counts - instead of a list of observed outcomes we have a vector where the ith entry is the number of times outcome i was observed
    //outcomes are assumed to be in lexicographic order
    vector<int> counts(pow(2,nvar),0);
    for(int i=0;i<D;i++){
        counts[data[i]]+=1;
    }
    //next we ensure that the number of comparisons is not more than the maximum possible
    long long int maxComp= (pow(2,nvar)*(pow(2,nvar)-1))/2;
    if(nTests>maxComp){
        nTests=maxComp;
    }
    //Comparisons will be a 0/1 vector recording which outcome pairs have been done already.
    //Every outcome corresponds to a number between 0 and 2^n-1 (they are enumerated lexicographically)
    //A distinct outcome pair can thus be written (a,b) a=/=b. If we order outcome pairs in lexicographic order (when considered as coordinate pair) then the ithe pair corresponds to ith entry of comparisons
    //That is (1,0), (2,0), (2,1),(3,0),(3,1),(3,2),(4,0),........
    vector<int> Comparisons(maxComp,0);
    //if nTests=maxposns then we can just move through the outcome pairs systematically as all need to be tested
    //Otherwise, we must select a (not yet considered) outcome pair at random each time
    //outposn1 - lexicographic position of outcome 1 in our pair
    //outposn2 - lexicographic position of outcome 2
    long long int outPosn1=1;
    long long int outPosn2=0;
    //LexMult is a vector of lexicographic multipliers that we use to move between an outcome and its lex position
    vector<long long int> LexMult(nvar);
    for(int i=0;i<nvar;i++){
        LexMult[i]=pow(2,nvar-i-1);
    }
    //initiate a random generator for uniform selection of an outcome pair 
    //randomly selects a number between 0 and  maxComp-1, this then gives an outcome pair by the enumeration of pairs we use in Comparisons
    struct timeval tv;
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
        typedef std::mt19937 G;
        G g(mySeed);
        typedef std::uniform_int_distribution<long long int> Dist;
        Dist uni(0, maxComp-1);
    //later in the function we need to perform dominance testing. We must add 1 to our breaks vector before feeding the CP-net into this function
    //This is an issue left over from translating functions from R to C++
    vector<int> breaks=cpnet.cptBreaks;
    for(int i=0;i<breaks.size();i++){
      breaks[i]+=1;
    }
    //supportedRelns will count the number of outcomes pairs we test where the cpnet does not contradict the data;
    long long int supportedRelns=0;
    //perform nTests many outcome pair tests:
    for(long long int i=0;i<nTests;i++){
        //if nTests<maxposns, we need to randomly select a new outcome pair.
        if(nTests<maxComp){
            bool newOutcome=false;
            while(!newOutcome){
                //randomly select an outcome pair by generating its lexicographic position
                long long int outPair=uni(g);
                //Check Comparisons vector to see if it has been considered previously
                if(Comparisons[outPair]==0){
                    //If it hasn't been considered previously, mark in Comparisons that we have now done this pair
                    Comparisons[outPair]=1;
                    newOutcome=true;
                    //convert the lexicographic position of the pair into the coordinate pair that gives the lex positions of outcome 1 and 2 in the pair
                    outPosn1=2;
                    bool identified=false;
                    while(!identified){
                        if(outPair<(outPosn1*(outPosn1-1))/2){
                            identified=true;
                            outPosn1-=1;
                        }
                        else{
                            outPosn1+=1;
                        }
                    }
                    outPosn2=outPair - ((outPosn1*(outPosn1-1))/2);
                }
                //If it has been considered previously, generate a new pair and try again
            }
        }
        //If we are systematically moving through the pairs, then mark this pair as considered
        if(nTests==maxComp){
            Comparisons[i]=1;
        }
        // let us copy outPosn1/2 to new ints to preserve them
        long long int posn1=outPosn1;
        long long int posn2=outPosn2;
        //next we need to turn the pair of outcome positions into two outcomes (vectors out1 and out2)
        vector<int> out1(nvar,1);
        vector<int> out2(nvar,1);
        for(int j=0;j<nvar;j++){
            int div1=posn1/LexMult[j];
            if(div1==1){
                out1[j]=2;
                posn1-=LexMult[j];
            }
            int div2=posn2/LexMult[j];
            if(div2==1){
                out2[j]=2;
                posn2-=LexMult[j];
            }
        }
        //next we want to find the number of observed data points for the selected outcomes:
        int data1=counts[outPosn1];
        int data2=counts[outPosn2];
        if(data1>data2){
            //out1 preferred to out2 is only contradicted by the cpnet if it entails out2>out1 so perform dominance test out2>out1
            List DTOut=RSDQRankPriority(cpnet.structure,cpnet.domains,cpnet.cptEntries,breaks,out2,out1);
            bool outcome=DTOut.result;
            if(!outcome){
                //If dominance test false then cpnet does not entail out2>out1 so the data is not contradicted
                //Thus, we say the CP-net is consistent with this relation in the data and we add a support count
                supportedRelns+=1;
            }
        }
        if(data2>data1){
            //out2 preferred to out1 is only contradicted by the cpnet if it entails out1>out2 so perform dominance test out1>out2
            List DTOut=RSDQRankPriority(cpnet.structure,cpnet.domains,cpnet.cptEntries,breaks,out1,out2);
            bool outcome=DTOut.result;
            if(!outcome){
                //If dominance test false then cpnet does not entail out1>out2 so the data is not contradicted
                //Thus, we say the CP-net is consistent with this relation in the data and we add a support count
                supportedRelns+=1;
            }
        }
        if(data1==data2){
            //out1 equally preferred to out2 is contradicted by the cpnet if it entails either out1>out2 or out2>out1 so  dominance test both
            //First we dominance test out1>out2
            List DTOut=RSDQRankPriority(cpnet.structure,cpnet.domains,cpnet.cptEntries,breaks,out1,out2);
            bool outcome=DTOut.result;
            if(!outcome){
                //If dominance test is false, we then test out2>out1
                DTOut=RSDQRankPriority(cpnet.structure,cpnet.domains,cpnet.cptEntries,breaks,out2,out1);
                outcome=DTOut.result;
                if(!outcome){
                    //If both dominance tests are false, neither direction is entailed and the CP-net is consistent with the data relation so we add a support count
                    supportedRelns+=1;
                }
            }
        }
        //if we are moving through the pairs systematically, we then move on to next pair:
        if(nTests==maxComp){
            if(outPosn2==(outPosn1-1)){
                outPosn1+=1;
                outPosn2=0;
            }
            else{
                outPosn2+=1;
            }
        }
    }
    //supportedRelns now counts the number of tests (outcome pairs) for which the CP-net was consistent with the data relation
    //DOC is then obtained by dividing by the number of tests performed
    return (double)supportedRelns/(double)nTests;
}

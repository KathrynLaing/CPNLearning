#include "10_DirichletSample.h"
#include "cpnInput.h"
#include "cpnOutput.h"
#include <math.h>
#include <algorithm> 
#include <random>
#include <ctime>
#include <chrono> 
#include <cmath>

LearningOutput cpnLearn(int nvar, std::vector<int> data, std::vector<double> beta, int DataSize, int mcReps, double signif, bool StartEmpty=true, vector<int> StartStructure={0}, string CpnSeqFilename="/file/path/output.csv"){
    // nvar - number of variables
    // data - vector of observed outcome choices (outcomes enumerated as integers indexed from 0)
    // beta - prior parameters
    // mcReps - Size of the Monte Carlo sample
    // signif - Change threshold in percent eg 5 -> 5%
    // StartEmpty - true if we want to use an empty as the starting structure
    // StartStructure - the starting structure if specified (i.e. StartEmpty=false), input as flattened adjacency matrix. Set to {0} if we are using an empty start.
    // CpnSeqFilename - file to write the sequence of CPNs visited by learning procedure to

    // Create initial CPN called cpnet that is the empty CPN if StartEmpty=true and the specified structure otherwise
    // Log the start time of the procedure
    auto start = std::chrono::high_resolution_clock::now();
    vector<double> StartBeta=beta;
    struct cpn cpnet;
    cpnet.size=nvar;
    std::vector<int> AdjMat(nvar*nvar,0);
    //if StartEmpty=false set structure as specified
    if(!StartEmpty){
        AdjMat=StartStructure;
    }
    cpnet.structure=AdjMat;
    // Set all domain sizes to 2
    std::vector<int> nDom(nvar,2);
    cpnet.domains=nDom;
    // Initialise the CPTs appropriate for the starting structure and fill each row with the rule "1 preferred to 2"
    std::vector<int> breaks(nvar+1);
    breaks[0]=0;
    if(StartEmpty){
        for(int i=0;i<nvar;i++){
            breaks[i+1]=breaks[i]+2;
        }
    }
    else{
        for(int i=0;i<nvar;i++){
        	// nPa= (#parents of variable i) +1
            int nPa=1;
            for(int j=0;j<nvar;j++){
                nPa+=AdjMat[j*nvar+i];
            }
            breaks[i+1]=breaks[i]+pow(2,nPa);
        }
    }
    cpnet.cptBreaks=breaks;
    std::vector<int> entries(breaks[nvar],1);
    for(int i=0; i<breaks[nvar]/2; i++){
        entries[2*i+1]=2;
    }
    cpnet.cptEntries=entries;
    //scores - the CPT scores of the variables
    std::vector<double> scores(nvar);
    //delta - matrix of delta values for all possible edge changes
    double delta[nvar][nvar];
    //cyclic - cyclicity matrix. 0 if an edge can be added, -1 if it cant (i.e. adding this edge creates cycles)
    int cyclic[nvar][nvar];
    // cyclic is a 0 matrix with -1 on diagonal if starting structure is empty
    for(int i=0;i<nvar;i++){
        for(int j=0;j<nvar;j++){
            if(i==j){
                cyclic[i][j]=-1;
            }
            else{
                cyclic[i][j]=0;
            }
        }
    }
    // Calculate the cyclicity matrix if starting structure is non empty
    if(!StartEmpty){
        for(int i=0;i<nvar;i++){
            //For each variable identify its ancestors
            vector <int> Anc(nvar,0);
            int check=0;
            int nAnc=0;
            for(int j=0;j<nvar;j++){
                Anc[j]=AdjMat[j*nvar+i];
                nAnc+=Anc[j];
            }
            while(nAnc!=check){
                check=nAnc;
                for(int j=0;j<nvar;j++){
                    if(Anc[j]==1){
                        for(int k=0;k<nvar;k++){
                            if(Anc[k]==0){
                                Anc[k]=AdjMat[k*nvar+j];
                                nAnc+=Anc[k];
                            }
                        }
                    }
                }
            }
            //Enter -1 to cyclic to show we cant add an edge from a variable to its ancestor as this would be a cycle
            //(note that such edges cannot already be present as the starting structure is acyclic and cannot be added)
            for(int j=0;j<nvar;j++){
                if(Anc[j]==1){
                    cyclic[i][j]=-1;
                }
            }
            //All other edges can be added/removed without introducing cycles
        }
    }
    // The prospective CPTs are a vector of vectors. The ith entry corresponds the edge A->B located in the ith position when adjacency matrices are written as vectors (as in our CP-net structures and input starting structures)
    // Equivalently, the ith entry is the A->B that comes ith lexicographically if edges are considered as coordinates (A,B)
    // This ith entry gives the optimal CPT(B) for the scenario where A->B is changed (add/removed) in our structure.
    // If the edge A->B is changed at some point, we use this entry to update our CPTs so that they remain optimal
    vector<vector<int> > futureCPTs(nvar*nvar);
    // nOut is the number of outcomes
    int nOut=beta.size();
    // LexMult is a lexicographic multiplier vector
    std::vector<int> LexMult(nvar);
    LexMult[nvar-1]=1;
    for(int i=0;i<nvar-1;i++){
        LexMult[(nvar-2-i)]=2*LexMult[(nvar-1-i)];
    }
    //LearnedCpns will contain the CP-nets obtained after each edge change
    vector<cpn> LearnedCpns;
    // update beta observed data points
    for(int j=0; j<DataSize; j++){
        beta[data[j]]+=1;
    }
    // Generate the Dirichlet sample:
    std::vector<double> FullSample;
    FullSample=DirSample(beta,nOut,mcReps);
    // for each variable calulate the optimal CPT and the associated table score
    for(int j=0; j<nvar; j++){
        // Let Pa be a 0/1 vector of parents, let nPa be the number of parents
        std::vector<int> Pa(nvar,0);
        int nPa=0;
        for(int k=0;k<nvar;k++){
            if(cpnet.structure[k*nvar+j]==1){
                Pa[k]=1;
                nPa+=1;
            }
        }
        // nOther is the number of variables not X or parents of X
        int nOther=nvar-1-nPa;
        // We use lexicographic ordering to encode each outcome as an integer.
        // We need to consider, for each Pa(X) assignment, every possible outcome in which X=x1 and every possible outcome in which X=x2
        // OtherWeights is a list of the "weights" each assignment to the other varaibles can add to the lexicographic value of an outcome
        // otherAssts cycles through the possible addignments to the other variables
        std::vector<int> otherWeights(pow(2,nOther));
        std::vector<int> otherAssts(nOther,0);
        std::vector<int> otherMult(nOther);
        if(nOther>0){
            int counter=0;
            for(int k=0;k<nvar;k++){
                if(k!=j){
                    if(Pa[k]==0){
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
        // for each dirichlet draw determine which CPT is supported
        // cptsSupported will store the UNIQUE supported CPTs
        // cptWeight will be the # of 1s in the CPTs in cptsSupported
        // cptsSupportCount will be the number of supporting draws for these CPTs
        std::vector<int> cptsSupported={};
        std::vector<int> cptWeight={};
        std::vector<int> cptsSupportCount={};
        // for each dirichlet sample:
        for(int k=0;k<mcReps;k++){
            // sample is the currently considered dirichlet draw
            std::vector<double> sample(nOut);
            for(int l=0;l<nOut;l++){
                sample[l]=FullSample[k*nOut + l];
            }
            // agreement gives, for each parental assignment, which preference rule sample supports - 0 for 1>2, 1 for 2>1, 2 for both (i.e. equal sums)
            std::vector<int> agreement(pow(2,nPa));
            // paAsst cycles through the possible parental assignments
            std::vector<int> paAsst(nPa,0);
            for(int l=0;l<pow(2,nPa);l++){
                // calculate the fixed weight added to lexicographic position by the parental assignment:
                int fixedWeight=0;
                int counter=0;
                for(int m=0;m<nvar;m++){
                    if(Pa[m]==1){
                        fixedWeight+=paAsst[counter]*LexMult[m];
                        counter+=1;
                    }
                }
                // group1 gives the lexicographic positions of the outcomes with this parent asst and variable j=0
                // group2 gives the ones with variable j=1
                std::vector<int> group1(pow(2,nOther));
                std::vector<int> group2(pow(2,nOther));
                if(nOther==0){
                    group1[0]=fixedWeight;
                    group2[0]=fixedWeight+LexMult[j];
                }
                else{
                    for(int m=0;m<pow(2,nOther);m++){
                        group1[m]=otherWeights[m]+fixedWeight;
                        group2[m]=otherWeights[m]+fixedWeight+LexMult[j];
                    }
                }
                // value1/2 gives the sum sample values for the outcomes of group1/2
                double value1=0.0;
                double value2=0.0;
                for(int m=0;m<pow(2,nOther);m++){
                    value1+=sample[group1[m]];
                    value2+=sample[group2[m]];
                }
                //now that we know the values we can determine which rule(s) the sample supports and record in the agreement vector
                if(value1>value2){
                    agreement[l]=0;
                }
                else{
                    if(value2>value1){
                        agreement[l]=1;
                    }
                    else{
                        agreement[l]=2;
                    }
                }
                //go to the next pa Asst
                if(l<pow(2,nPa)-1){
                    int max=0;
                    for(int m=0;m<nPa;m++){
                        if(paAsst[m]==0){
                            max=m;
                        }
                    }
                    paAsst[max]=1;
                    for(int m=max+1;m<nPa;m++){
                        paAsst[m]=0;
                    }
                }
            }
            //agreement now gives the cpt supported by sample
            //however, if agreement has any 2 entries, then more than one CPT is supported and we must randomly select one
            //we randomly select one by replacing each case of a 2 entry by 0 or 1 randomly (if both rule directions equally supported, randomly assign a rule)
            //First, create a random generator that assigns the values of 0 or 1 with equal probability
            struct timeval tv;
            gettimeofday(&tv,0);
            unsigned long mySeed = tv.tv_sec + tv.tv_usec;
            std::default_random_engine eng(mySeed);
            std::uniform_int_distribution<> UnifIntDist(0, 1);
            //Now replace each 2 in agreement with a 0 or 1 at random
            for(int l=0;l<pow(2,nPa);l++){
                if(agreement[l]==2){
                    agreement[l]=UnifIntDist(eng);
                }
            }
            // We now have the single CPT supported by the dirichlet sample draw
            // let us add this to the list of supported CPTs if its new or add a support count if not
            // first determine its weight (number of 1s in this CPT)
            int weight=0;
            for(int l=0;l<pow(2,nPa);l++){
                weight+=agreement[l];
            }
            // how many CPTs are supported so far?
            int nSupCpts=cptWeight.size();
            if(nSupCpts==0){
                // if no supported CPTs yet, add this one, its weight and 1 support count to the vectors respectively
                for(int l=0;l<pow(2,nPa);l++){
                    cptsSupported.push_back(agreement[l]);
                }
                cptWeight.push_back(weight);
                cptsSupportCount.push_back(1);
            }
            else{
            	// If there are some supported CPTs already:
                //let us see if our CPT matches any of the others
                bool NewCPT=true;
                for(int l=0;l<nSupCpts;l++){
                    //for each existing supported CPT, see if it has same weight
                    if(cptWeight[l]==weight){
                        //if they have equal weight, see if they're the same
                        //ndiffs - number of differences in their entries
                        int ndiffs=0;
                        for(int m=0;m<pow(2,nPa);m++){
                            if(agreement[m]!=cptsSupported[(pow(2,nPa)*l)+m]){
                               ndiffs+=1;
                               break; 
                            }
                        }
                        if(ndiffs==0){
                            //ndiffs=0 means that our CPT matches one already in the list, so we add a support count and stop looking
                            cptsSupportCount[l]+=1;
                            NewCPT=false;
                            break;
                        }
                    }
                }
                if(NewCPT){
                    //if we didnt find it in the list, add it to the end, add its weight, and give it a support count of 1
                    for(int l=0;l<pow(2,nPa);l++){
                        cptsSupported.push_back(agreement[l]);
                    }
                    cptWeight.push_back(weight);
                    cptsSupportCount.push_back(1);
                }
            }
        }
        //cptsSupported now gives the unique CPTs supported by the Dirichlet samples
        //cptsSupportCount gives their support counts, which is maximal?
        // max - maximal support count
        // maxposn - the index of the maximally supported CPT in the list of unique CPTs
        int max=0;
        int maxposn;
        int nSupCpts=cptsSupportCount.size();
        for(int k=0; k<nSupCpts; k++){
            if(cptsSupportCount[k]>max){
                max=cptsSupportCount[k];
                maxposn=k;
            }
        }
        //Let us identify ALL SupportedCPTs with this score and randomly select one to assign to the CPN (there may be multiple maximally supported CPTs)
        // nBestCPTs - number of maximally supported CPTs
        // BestCPTs - the positions of these CPTs within the cptsSupported vector
        int nBestCPTs=1;
        vector<int> BestCPTs={maxposn};
        for(int k=maxposn+1; k<nSupCpts; k++){
            if(cptsSupportCount[k]==max){
                BestCPTs.push_back(k);
                nBestCPTs+=1;
            }
        }
        // Random generator to randomly sample a number between 0 and nBestCPTs-1
        struct timeval tv;
        gettimeofday(&tv,0);
        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
        std::default_random_engine eng(mySeed);
        std::uniform_int_distribution<> UnifIntDist(0, nBestCPTs-1);
        //Pick the "optimal" CPT by randomly selecting one of the maximally supported
        int ChosenCPTPosn=BestCPTs[UnifIntDist(eng)];
        // The scorre of  Variable j is then the number of support counts/ size of the dirichlet sample
        scores[j]=(double) cptsSupportCount[ChosenCPTPosn]/(double) mcReps;
        // assign variable j the optimal CPT by editing the variable j section of cptEntries of our cp-net
        breaks=cpnet.cptBreaks;
        std::vector<int> cpt=cpnet.cptEntries;
        int counter=0;
        for(int k=0; k<pow(2,nPa); k++){
        	// A 0 entry in optimal CPT implies the preference rule 1>2
            if(cptsSupported[ChosenCPTPosn*pow(2,nPa)+k] == 0){
                cpt[breaks[j]+counter]=1;
                cpt[breaks[j]+counter+1]=2;
            }
            else{// A 1 entry in optimal CPT implies the preference rule 2>1
                cpt[breaks[j]+counter]=2;
                cpt[breaks[j]+counter+1]=1;
            }
            counter+=2;
        }
        cpnet.cptEntries=cpt;
    }
    // We have now calculated the scores of all variables and assigned all variables their optimal CPTs
    // Now we must calculate the delta values
    // Simultaneously we identify and record the potential new CPTs in the futureCPTs matrix (these are the CPT(Y)s that would be optimal if the edge X->Y is changed)
    AdjMat=cpnet.structure;
    //For each edge we can add/remove, calculate the delta value
    for(int j=0;j<nvar;j++){
        for(int k=0;k<nvar;k++){
        	// we are here considering the edge from variable j -> variable k
            // X->X edges can be ignored as we are only ocnsidering acyclic structures
            if(j!=k){
                // If it is an edge for removal (already in the structure)
                if(AdjMat[j*nvar+k]==1){
                    // Let Pa be a 0/1 vector of the parents of k EXCLUDING j
                    // nPa be the number of parents in Pa
                    std::vector<int> Pa(nvar,0);
                    int nPa=0;
                    for(int m=0;m<nvar;m++){
                        if(m!=j){
                            if(cpnet.structure[m*nvar+k]==1){
                                Pa[m]=1;
                                nPa+=1;
                            }
                        }
                    }
                    // We now proceed to calculate the score for k if it had parents Pa (i.e if the j->k edge had been removed) in the same way as above.
                    // This will also generate an optimal CPT for this scenario.
                    // nother is the number of variables that are not k or in Pa
                    int nOther=nvar-1-nPa;
                    // otherWeights - list of all possible weights that the other variables can add to lexicographic posn
                    std::vector<int> otherWeights(pow(2,nOther));
                    std::vector<int> otherAssts(nOther,0);
                    std::vector<int> otherMult(nOther);
                    if(nOther>0){
                        int counter=0;
                        for(int m=0;m<nvar;m++){
                            if(m!=k){
                                if(Pa[m]==0){
                                    otherMult[counter]=LexMult[m];
                                    counter+=1;
                                }
                            }
                        }
                        for(int m=0;m<pow(2,nOther);m++){
                            int weight=0;
                            for(int n=0;n<nOther;n++){
                                weight+= otherAssts[n]*otherMult[n];
                            }
                            otherWeights[m]=weight;
                            if(m!=(pow(2,nOther)-1)){
                                int max=0;
                                for(int n=0;n<nOther;n++){
                                    if(otherAssts[n]==0){
                                        max=n;
                                    }
                                }
                                otherAssts[max]=1;
                                for(int n=max+1;n<nOther;n++){
                                    otherAssts[n]=0;
                                }
                            }
                        }
                    }
                    // for each sampled draw, let us record the supported CPT
                    // cptsSupported will store the UNIQUE supported CPTs
                    // cptWeight will be their weights (# of 1s)
                    // cptsSupportCount will be the number of supporting draws 
                    std::vector<int> cptsSupported={};
                    std::vector<int> cptWeight={};
                    std::vector<int> cptsSupportCount={};
                    // for each sample:
                    for(int l=0;l<mcReps;l++){
                        // sample is the currently considered dirichlet draw
                        std::vector<double> sample(nOut);
                        for(int m=0;m<nOut;m++){
                            sample[m]=FullSample[l*nOut + m];
                        }
                        // agreement gives, for each assignment to tha Pa variables, which rule sample supports - 0 for 1>2, 1 for 2>1, 2 for both
                        std::vector<int> agreement(pow(2,nPa));
                        // paAsst cycles through the assignments to Pa (the parents of k excluding j)
                        std::vector<int> paAsst(nPa,0);
                        for(int m=0;m<pow(2,nPa);m++){
                            // calculate the fixed value added to position number by the parental assignment:
                            int fixedWeight=0;
                            int counter=0;
                            for(int n=0;n<nvar;n++){
                                if(Pa[n]==1){
                                    fixedWeight+=paAsst[counter]*LexMult[n];
                                    counter+=1;
                                }
                            }
                            // group1 gives the positions of the outcomes with this parent asst and var k=0
                            // group2 gives the ones with var k=1
                            std::vector<int> group1(pow(2,nOther));
                            std::vector<int> group2(pow(2,nOther));
                            if(nOther==0){
                                group1[0]=fixedWeight;
                                group2[0]=fixedWeight+LexMult[k];
                            }
                            else{
                                for(int n=0;n<pow(2,nOther);n++){
                                    group1[n]=otherWeights[n]+fixedWeight;
                                    group2[n]=otherWeights[n]+fixedWeight+LexMult[k];
                                }
                            }
                            // value1/2 give the sum of the values of group1/2 outcomes in the sample
                            double value1=0.0;
                            double value2=0.0;
                            for(int n=0;n<pow(2,nOther);n++){
                                value1+=sample[group1[n]];
                                value2+=sample[group2[n]];
                            }
                            //now that we know the values we can put which rule(s) the sample supports into the agreement vector
                            if(value1>value2){
                                agreement[m]=0;
                            }
                            else{
                                if(value2>value1){
                                    agreement[m]=1;
                                }
                                else{
                                    agreement[m]=2;
                                }
                            }
                            //go to the next pa Asst
                            if(m<pow(2,nPa)-1){
                                int max=0;
                                for(int n=0;n<nPa;n++){
                                    if(paAsst[n]==0){
                                        max=n;
                                    }
                                }
                                paAsst[max]=1;
                                for(int n=max+1;n<nPa;n++){
                                    paAsst[n]=0;
                                }
                            }
                        }
                        //agreement now encodes which cpts are supported by sample
                        //let us now randomly select one by replacing each case of 2 by 0 or 1 randomly
                        struct timeval tv;
                        gettimeofday(&tv,0);
                        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                        std::default_random_engine eng(mySeed);
                        std::uniform_int_distribution<> UnifIntDist(0, 1);
                        for(int m=0;m<pow(2,nPa);m++){
                            if(agreement[m]==2){
                                agreement[m]=UnifIntDist(eng);
                            }
                        }
                        // let us add this to the list of supported CPTs if its new or add a support count if not
                        // first determine its weight
                        int weight=0;
                        for(int m=0;m<pow(2,nPa);m++){
                            weight+=agreement[m];
                        }
                        // how many CPTs are supported so far?
                        int nSupCpts=cptWeight.size();
                        if(nSupCpts==0){
                            // if no supported CPTs yet, add this one, its weight and 1 support count to the vectors respectively
                            for(int m=0;m<pow(2,nPa);m++){
                                cptsSupported.push_back(agreement[m]);
                            }
                            cptWeight.push_back(weight);
                            cptsSupportCount.push_back(1);
                        }
                        else{
                            //let us see if our CPT matches any of the others
                            bool NewCPT=true;
                            for(int m=0;m<nSupCpts;m++){
                                //for each existing supported CPT, see if it has same weight
                                if(cptWeight[m]==weight){
                                    //if they have equal weight, see if they're the same
                                    int ndiffs=0;
                                    for(int n=0;n<pow(2,nPa);n++){
                                        if(agreement[n]!=cptsSupported[(pow(2,nPa)*m)+n]){
                                           ndiffs+=1;
                                           break; 
                                        }
                                    }
                                    if(ndiffs==0){
                                        //if we have found the CPT in the list then add a support count and stop looking
                                        cptsSupportCount[m]+=1;
                                        NewCPT=false;
                                        break;
                                    }
                                }
                            }
                            if(NewCPT){
                                //if we didnt find it in the list, add it to the end, add its weight, add a support count of 1
                                for(int m=0;m<pow(2,nPa);m++){
                                    cptsSupported.push_back(agreement[m]);
                                }
                                cptWeight.push_back(weight);
                                cptsSupportCount.push_back(1);
                            }
                        }
                    }
                    //cptsSupported now gives the CPTs supported by the Dirichlet sample
                    //cptsSupportCount gives the support counts, which is max?
                    int max=0;
                    int maxposn;
                    int nSupCpts=cptsSupportCount.size();
                    for(int l=0; l<nSupCpts; l++){
                        if(cptsSupportCount[l]>max){
                            max=cptsSupportCount[l];
                            maxposn=l;
                        }
                    }
                    //Let us identify ALL SupportedCPTs with maximal score and randomly select one to assign as the (potential) optimal CPT
                    int nBestCPTs=1;
                    vector<int> BestCPTs={maxposn};
                    for(int l=maxposn+1; l<nSupCpts; l++){
                        if(cptsSupportCount[l]==max){
                            BestCPTs.push_back(l);
                            nBestCPTs+=1;
                        }
                    }
                    struct timeval tv;
                    gettimeofday(&tv,0);
                    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                    std::default_random_engine eng(mySeed);
                    std::uniform_int_distribution<> UnifIntDist(0, nBestCPTs-1);
                    int ChosenCPTPosn=BestCPTs[UnifIntDist(eng)];
                    // Delta value of the edge j->k = table score of variable k with this new optimal CPT/current score for variable k
                    delta[j][k]=(double) cptsSupportCount[ChosenCPTPosn]/((double) mcReps * scores[k]);
                    // identify the (randomly selected) optimal CPT and assign to future CPTs
                    breaks=cpnet.cptBreaks;
                    std::vector<int> cpt(pow(2,nPa+1));
                    int counter=0;
                    for(int l=0; l<pow(2,nPa); l++){
                        if(cptsSupported[(ChosenCPTPosn*pow(2,nPa))+l] == 0){
                            cpt[counter]=1;
                            cpt[counter+1]=2;
                        }
                        else{
                            cpt[counter]=2;
                            cpt[counter+1]=1;
                        }
                        counter+=2;
                    }
                    futureCPTs[j*nvar+k].resize(pow(2,nPa+1));
                    futureCPTs[j*nvar+k]=cpt;
                }
                else{// if it is an edge to be added (i.e. j->k is not already in the structure)
                    // Let Pa be a 0/1 vector of the parents of k plus j
                    // nPa - number of parents in Pa
                    std::vector<int> Pa(nvar,0);
                    int nPa=0;
                    for(int m=0;m<nvar;m++){
                        if(m!=j){
                            if(cpnet.structure[m*nvar+k]==1){
                                Pa[m]=1;
                                nPa+=1;
                            }
                        }
                        else{
                            Pa[m]=1;
                            nPa+=1;
                        }
                    }
                    // We now calculate the score for variable k if it had parents Pa (i.e. the score of k if edge j->k was added)
                    // This will also determine the optimal CPT in this scenario
                    // We calculate this in the same manner as above
                    // nother is the number of variables not k or parents of k or j
                    int nOther=nvar-1-nPa;
                    // otherWeights - list of all possible weights the other variables can add to lexicographic posn
                    std::vector<int> otherWeights(pow(2,nOther));
                    std::vector<int> otherAssts(nOther,0);
                    std::vector<int> otherMult(nOther);
                    if(nOther>0){
                        int counter=0;
                        for(int m=0;m<nvar;m++){
                            if(m!=k){
                                if(Pa[m]==0){
                                    otherMult[counter]=LexMult[m];
                                    counter+=1;
                                }
                            }
                        }
                        for(int m=0;m<pow(2,nOther);m++){
                            int weight=0;
                            for(int n=0;n<nOther;n++){
                                weight+= otherAssts[n]*otherMult[n];
                            }
                            otherWeights[m]=weight;
                            if(m!=(pow(2,nOther)-1)){
                                int max=0;
                                for(int n=0;n<nOther;n++){
                                    if(otherAssts[n]==0){
                                        max=n;
                                    }
                                }
                                otherAssts[max]=1;
                                for(int n=max+1;n<nOther;n++){
                                    otherAssts[n]=0;
                                }
                            }
                        }
                    }
                    // for each sampled draw, let us record the supported CPT
                    // cptsSupported will store the UNIQUE supported CPTs
                    // cptWeight will be their weights (# of 1s)
                    // cptsSupportCount will be the number of supporting draws 
                    std::vector<int> cptsSupported={};
                    std::vector<int> cptWeight={};
                    std::vector<int> cptsSupportCount={};
                    // for each sample:
                    for(int l=0;l<mcReps;l++){
                        // sample is the currently considered sample value set
                        std::vector<double> sample(nOut);
                        for(int m=0;m<nOut;m++){
                            sample[m]=FullSample[l*nOut + m];
                        }
                        // agreement gives, for each assignment to Pa, which rule sample supports - 0 for 1>2, 1 for 2>1, 2 for both
                        std::vector<int> agreement(pow(2,nPa));
                        // paAsst cycles through all possible Pa assignments
                        std::vector<int> paAsst(nPa,0);
                        for(int m=0;m<pow(2,nPa);m++){
                            // calculate the fixed value added to lexicographic position number by the parental assignment:
                            int fixedWeight=0;
                            int counter=0;
                            for(int n=0;n<nvar;n++){
                                if(Pa[n]==1){
                                    fixedWeight+=paAsst[counter]*LexMult[n];
                                    counter+=1;
                                }
                            }
                            // group1 gives the lexicographic positions of the outcomes with this parent asst and var k=0
                            // group2 gives the ones with var k=1
                            std::vector<int> group1(pow(2,nOther));
                            std::vector<int> group2(pow(2,nOther));
                            if(nOther==0){
                                group1[0]=fixedWeight;
                                group2[0]=fixedWeight+LexMult[k];
                            }
                            else{
                                for(int n=0;n<pow(2,nOther);n++){
                                    group1[n]=otherWeights[n]+fixedWeight;
                                    group2[n]=otherWeights[n]+fixedWeight+LexMult[k];
                                }
                            }
                            // value1/2 give the sum of the values of the outcomes in group1/2 in the sample
                            double value1=0.0;
                            double value2=0.0;
                            for(int n=0;n<pow(2,nOther);n++){
                                value1+=sample[group1[n]];
                                value2+=sample[group2[n]];
                            }
                            //now that we know the sum values we can put which rule(s) the sample supports into the agreement vector
                            if(value1>value2){
                                agreement[m]=0;
                            }
                            else{
                                if(value2>value1){
                                    agreement[m]=1;
                                }
                                else{
                                    agreement[m]=2;
                                }
                            }
                            //go to the next pa Asst
                            if(m<pow(2,nPa)-1){
                                int max=0;
                                for(int n=0;n<nPa;n++){
                                    if(paAsst[n]==0){
                                        max=n;
                                    }
                                }
                                paAsst[max]=1;
                                for(int n=max+1;n<nPa;n++){
                                    paAsst[n]=0;
                                }
                            }
                        }
                        //agreement now encodes WHICH cpts are supported by sample
                        //let us now randomly select one by replacing each 2 entry by 0 or 1 randomly
                        struct timeval tv;
                        gettimeofday(&tv,0);
                        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                        std::default_random_engine eng(mySeed);
                        std::uniform_int_distribution<> UnifIntDist(0, 1);
                        for(int m=0;m<pow(2,nPa);m++){
                            if(agreement[m]==2){
                                agreement[m]=UnifIntDist(eng);
                            }
                        }
                        // let us add this to the list of supported CPTs if its new or add a support count if not
                        // first determine its weight (# of 1s)
                        int weight=0;
                        for(int m=0;m<pow(2,nPa);m++){
                            weight+=agreement[m];
                        }
                        // how many CPTs are supported so far?
                        int nSupCpts=cptWeight.size();
                        if(nSupCpts==0){
                            // if no supported CPTs yet, add this one, its weight and 1 support count to the vectors respectively
                            for(int m=0;m<pow(2,nPa);m++){
                                cptsSupported.push_back(agreement[m]);
                            }
                            cptWeight.push_back(weight);
                            cptsSupportCount.push_back(1);
                        }
                        else{// if there are some supported CPTs already
                            //let us see if our CPT matches any of the others
                            bool NewCPT=true;
                            for(int m=0;m<nSupCpts;m++){
                                //for each existing supported CPT, see if it has same weight
                                if(cptWeight[m]==weight){
                                    //if they have equal weight, see if they're the same
                                    int ndiffs=0;
                                    for(int n=0;n<pow(2,nPa);n++){
                                        if(agreement[n]!=cptsSupported[(pow(2,nPa)*m)+n]){
                                           ndiffs+=1;
                                           break; 
                                        }
                                    }
                                    if(ndiffs==0){
                                        //if we have found it in the list then add a support count and stop looking
                                        cptsSupportCount[m]+=1;
                                        NewCPT=false;
                                        break;
                                    }
                                }
                            }
                            if(NewCPT){
                                //if we didnt find it in the list, add it to the end, add the CPT weight, and a support count of 1
                                for(int m=0;m<pow(2,nPa);m++){
                                    cptsSupported.push_back(agreement[m]);
                                }
                                cptWeight.push_back(weight);
                                cptsSupportCount.push_back(1);
                            }
                        }
                    }
                    //cptsSupported now gives the CPTs supported by the Dirichlet sample
                    //cptsSupportCount gives the support counts of each, which is max?
                    int max=0;
                    int maxposn;
                    int nSupCpts=cptsSupportCount.size();
                    for(int l=0; l<nSupCpts; l++){
                        if(cptsSupportCount[l]>max){
                            max=cptsSupportCount[l];
                            maxposn=l;
                        }
                    }
                    //Let us identify ALL SupportedCPTs with maximal score and randomly select one to assign as the (potential) optimal CPT
                    int nBestCPTs=1;
                    vector<int> BestCPTs={maxposn};
                    for(int l=maxposn+1; l<nSupCpts; l++){
                        if(cptsSupportCount[l]==max){
                            BestCPTs.push_back(l);
                            nBestCPTs+=1;
                        }
                    }
                    struct timeval tv;
                    gettimeofday(&tv,0);
                    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                    std::default_random_engine eng(mySeed);
                    std::uniform_int_distribution<> UnifIntDist(0, nBestCPTs-1);
                    // The chosen optimal CPT is the best scoring CPT in the case where k has its current parents plus j (i.e. when edge j->k has been added)
                    int ChosenCPTPosn=BestCPTs[UnifIntDist(eng)];
                    // Delta value for j->k is table score for this optimal CPT/current variable k score
                    delta[j][k]=(double) cptsSupportCount[ChosenCPTPosn]/((double) mcReps * scores[k]);
                    // identify the (randomly selected) optimal CPT and assign to future CPTs - this is the optimal CPT(k) in the case where edge j->k is changed
                    breaks=cpnet.cptBreaks;
                    std::vector<int> cpt(pow(2,nPa+1));
                    int counter=0;
                    for(int l=0; l<pow(2,nPa); l++){
                        if(cptsSupported[(ChosenCPTPosn*pow(2,nPa))+l] == 0){
                            cpt[counter]=1;
                            cpt[counter+1]=2;
                        }
                        else{
                            cpt[counter]=2;
                            cpt[counter+1]=1;
                        }
                        counter+=2;
                    }
                    futureCPTs[j*nvar+k].resize(pow(2,nPa+1));
                    futureCPTs[j*nvar+k]=cpt;
                }
            }
        }
    }
    //We have now calculated the delta values for all possible edge changes. We have also recorded the optimal CPTs for all possible edge changes
    //Now that we have observed data and calculated all scores, deltas, optimal cpts, and prospective cpts we can commence learning
    //First, are there any valid edge changes we can make?
    //Continue will determine whether we can continue learning (i.e whether there are viable further edge changes to make)
    bool Continue=false;
    AdjMat=cpnet.structure;
    //Valid changes gives the delta value of valid changes and -1 for invalid changes
    std::vector<double> ValidChanges(nvar*nvar);
    //loop through all possible edges j->k
    for(int j=0;j<nvar;j++){
        for(int k=0;k<nvar;k++){
        	// If cyclic has entry -1 then the edge cannot be changed without introducing cycles so the change is invalid
            if(cyclic[j][k]!=0){
                ValidChanges[j*nvar +k]=-1;
            }
            else{
            	//depending on whether the edge is to be removed or added (in/not in the structure already), determine whether the delta is above the appropriate threshold
            	//if so, it is a valid edge change, enter delta value to ValidChanges. Otherwise, it is not a valid change, assign a -1 entry to ValidChanges
                if(AdjMat[j*nvar +k]==0){
                    if(delta[j][k]>(1+signif/100)){
                        ValidChanges[j*nvar +k]=delta[j][k];
                        Continue=true;
                    }
                    else{
                        ValidChanges[j*nvar +k]=-1;
                    }
                }
                else{
                    if(delta[j][k]>(1/(1+signif/100))){
                        ValidChanges[j*nvar +k]=delta[j][k];
                        Continue=true;
                    }
                    else{
                        ValidChanges[j*nvar +k]=-1;
                    }
                }
            }
        }
    }
    // If Continue=true, there must be at leat 1 valid change
    // While there are valid edge changes, make the best (largest delta) valid change and then update the structure, scores, and deltas and repeat
    // breaks and entries - the current CPTs
    breaks=cpnet.cptBreaks;
    entries=cpnet.cptEntries;
    while(Continue){
        // While there are valid edge changes to be made:
        // MaxDelta maximum delta of the valid changes
        // EdgeChange - index of (one of) the edge(s) with maximum delta
        // Edges are a number between 0 and nvar*nvar -1, it gives the position in the CP-net structure vector
        double MaxDelta=-1.0;
        int EdgeChange;
        for(int j=0; j<(nvar*nvar);j++){
            if(ValidChanges[j]>MaxDelta){
                MaxDelta=ValidChanges[j];
                EdgeChange=j;
            }
        }
        //Let us identify all valid edge changes with maximum delta (Best Edges) and select one at random to implement (ChosenEdgeChange):
        //As we are assessing equality of doubles we allow some leeway with respect to what constitutes equality (to allow for rounding error)
        int nBestEdges=1;
        vector<int> BestEdges={EdgeChange};
        for(int j=0; j<(nvar*nvar);j++){
            if(j!=EdgeChange){
                if(abs(MaxDelta-ValidChanges[j])<(double)0.5/(double)(mcReps*(mcReps-1))){
                    BestEdges.push_back(j);
                    nBestEdges+=1;
                }
            }
        }
        struct timeval tv;
        gettimeofday(&tv,0);
        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
        std::default_random_engine eng(mySeed);
        std::uniform_int_distribution<> UnifIntDist(0, nBestEdges-1);
        int ChosenEdgeChange=BestEdges[UnifIntDist(eng)];
        // If our chosen edge is X->Y then row gives X and col gives Y
        int row= ChosenEdgeChange/nvar;
        int col=ChosenEdgeChange-(row*nvar);
        //Implement this chosen edge change
        //Start by changing the structure as appropriate
        //Let remove tell us if it was an addition or a removal
        bool remove=false;
        AdjMat=cpnet.structure;
        //If X->Y not yet in structure, add the edge 
        if(AdjMat[ChosenEdgeChange]==0){
            AdjMat[ChosenEdgeChange]=1;
        }
        else{//If X->Y already present, remove it
            AdjMat[ChosenEdgeChange]=0;
            remove=true;
        }
        //update the structure in our CP-net
        cpnet.structure=AdjMat;
        //Update the score of variable Y by multiplying by the delta of X->Y
        scores[col]*=ValidChanges[ChosenEdgeChange];
        //update the CPT(Y) by changing it to the one stored in futureCPTs corresponding to this edge change
        //Start by updating the breaks vector
        std::vector<int> OldBreaks=cpnet.cptBreaks;
        std::vector<int> NewBreaks(nvar+1);
        //OldCptSize - size of CPT(Y) in previous CP-net
        int OldCptSize=OldBreaks[col+1]-OldBreaks[col];
        for(int j=0;j<nvar+1;j++){
            if(j<=col){
                NewBreaks[j]=OldBreaks[j];
            }
            else{
                if(remove){
                	//If we removed an edge, then Y lost a parent and the CPT size halved so the break points after CPT(Y) will all be smaller
                    NewBreaks[j]=OldBreaks[j]-(OldCptSize/2);
                }
                else{
                	//If we added an edge, then Y gained a parent and the CPT size doubled so the break points after CPT(Y) will all be larger
                    NewBreaks[j]=OldBreaks[j]+(OldCptSize);
                }
            }
        }
        cpnet.cptBreaks=NewBreaks;
        //Now update the entries
        //newcpt is the new CPT(Y) we want to insert
        vector<int> newcpt=futureCPTs[ChosenEdgeChange];
        //OldEntries - the entries in the current CP-net
        std::vector<int> OldEntries=cpnet.cptEntries;
        // NewEntreis - updated version
        std::vector<int> NewEntries(NewBreaks[nvar]);
        // Before CPT(Y), the entries are the same
        for(int j=0;j<OldBreaks[col];j++){
            NewEntries[j]=OldEntries[j];
        }
        if(remove){
        	//insert the new CPT(Y)
            for(int j=0;j<(OldCptSize/2);j++){
                NewEntries[OldBreaks[col]+j]=newcpt[j];
            }
            //After CPT(Y) it is the same entries as in OldEntries just shifted position due to the new CPT(Y) size [half previous size]
            for(int j=NewBreaks[col+1];j<NewBreaks[nvar];j++){
                NewEntries[j]=OldEntries[j+(OldCptSize/2)];
            }
        }
        else{
        	//insert the new CPT(Y)
            for(int j=0;j<(2*OldCptSize);j++){
                NewEntries[OldBreaks[col]+j]=newcpt[j];
            }
            //After CPT(Y) it is the same entries as in OldEntries just shifted position due to the new CPT(Y) size [double previous size]
            for(int j=NewBreaks[col+1];j<NewBreaks[nvar];j++){
                NewEntries[j]=OldEntries[j-(OldCptSize)];
            }
        }
        //replace entries vector in our CP-net
        cpnet.cptEntries=NewEntries;
        //The CP-net is now fully updated
        //LearnedCpns records each new CP-net we obtain after each edge change in order to see the behaviour of our algorithm.
        LearnedCpns.push_back (cpnet);
        //next we update the cyclicity matrix
        //If we add X->Y then adding Y->X is not possible as it would introduce cycles
        if(!remove){
            cyclic[col][row]=-1;
        }
        //anc is a 0/1 vector that has 1s for all ancestors of X and X itself
        vector<int> anc(nvar,0);
        int check=0;
        int nanc=0;
        //Add all parents of X
        for(int j=0;j<nvar;j++){
            anc[j]=AdjMat[j*nvar+row];
            nanc+=anc[j];
        }
        while(check!=nanc){
        	//Add all parents of variables in anc until there is no change
            check=nanc;
            for(int j=0;j<nvar;j++){
                if(anc[j]==1){
                    for(int k=0;k<nvar;k++){
                        if(anc[k]==0){
                            anc[k]=AdjMat[k*nvar+j];
                            nanc+=anc[k];
                        }
                    }
                }
            }
        }
        // Add X to anc
        anc[row]=1;
        //dec - descendants of col (Y) and Y (0/1 vector)
        //Calculated similarly to anc
        vector<int> dec(nvar,0);
        check=0;
        int ndec=0;
        for(int j=0;j<nvar;j++){
            dec[j]=AdjMat[col*nvar+j];
            ndec+=dec[j];
        }
        while(check!=ndec){
            check=ndec;
            for(int j=0;j<nvar;j++){
                if(dec[j]==1){
                    for(int k=0;k<nvar;k++){
                        if(dec[k]==0){
                            dec[k]=AdjMat[j*nvar+k];
                            ndec+=dec[k];
                        }
                    }
                }
            }
        }
        dec[col]=1;
        //For each edge k->j
        for(int j=0;j<nvar;j++){
            // calculate the descendants of j (dec2)
            check=0;
            ndec=0;
            vector<int> dec2(nvar,0);
            for(int k=0;k<nvar;k++){
                dec2[k]=AdjMat[j*nvar+k];
                ndec+=dec2[k];
            }
            while(check!=ndec){
                check=ndec;
                for(int k=0;k<nvar;k++){
                    if(dec2[k]==1){
                        for(int l=0;l<nvar;l++){
                            if(dec2[l]==0){
                                dec2[l]=AdjMat[k*nvar+l];
                                ndec+=dec2[l];
                            }
                        }
                    }
                }
            }
            for(int k=0;k<nvar;k++){
            	//Now consider each k->j edge
            	// k->k does not need updating as it is always impossible
                if(k!=j){
                    if(AdjMat[j*nvar+k]+AdjMat[k*nvar+j]==0){
                        //If k=/=j and they are unconnected (no edge either direction)
                        if(!remove){
                        	//If we have added an edge to the structure
                            if(cyclic[k][j]==0){
                            	//It is possible that edges that were previously possible to add may now introduce cycles.
                                //This is only true for k->j if j is an ancestor of X (or X) and k is a decsendant of Y (or Y)
                                if(anc[j]+dec[k]==2){
                                	//in this case, change cyclicity matrix to -1 to show that this edge is no longer possible to add
                                    cyclic[k][j]=-1;
                                }
                            }
                        }
                        else{
                        	//If we have removed an edge from the structure
                            if(cyclic[k][j]==-1){
                            	//It is possible that edges that were previously impossible to add may now be possible as they no longer produce cycles.
                                //This is only true for k->j if k is not a descendant of j in this new structure
                                if(dec2[k]==0){
                                	//in this case, change cyclicity matrix to 0 to show that this edge is now possible to add
                                    cyclic[k][j]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
        //cyclic is now an accurate cyclicity matrix for the updated CP-net structure
        //update the deltas (and associated prospective CPTs) for each edge of the form *->Y
        //For each variable j:
        for(int j=0;j<nvar;j++){
        	//Don't neet to consider Y->Y edge as this can't be changed (added)
            if(j!=col){
            	//We now update delta (and CPT) for the edge j->Y
                // If it is an edge for removal:
                if(AdjMat[j*nvar+col]==1){
                    // Let Pa be a 0/1 vector of parents of Y EXCLUDING j
                    // nPa - number of parents in Pa
                    std::vector<int> Pa(nvar,0);
                    int nPa=0;
                    for(int m=0;m<nvar;m++){
                        if(m!=j){
                            if(cpnet.structure[m*nvar+col]==1){
                                Pa[m]=1;
                                nPa+=1;
                            }
                        }
                    }
                    // We now calculate the score for Y (and optimal CPT(Y)) in the scenario where parent j is removed (i.e. when Y has parents Pa)
                    // We do this in the same way as score calculations were performed above
                    // nother is the number of variables not Y or in Pa
                    int nOther=nvar-1-nPa;
                    // otherWeights - list of all possible weights these other variables can add to lexicographic posn
                    std::vector<int> otherWeights(pow(2,nOther));
                    std::vector<int> otherAssts(nOther,0);
                    std::vector<int> otherMult(nOther);
                    if(nOther>0){
                        int counter=0;
                        for(int m=0;m<nvar;m++){
                            if(m!=col){
                                if(Pa[m]==0){
                                    otherMult[counter]=LexMult[m];
                                    counter+=1;
                                }
                            }
                        }
                        for(int m=0;m<pow(2,nOther);m++){
                            int weight=0;
                            for(int n=0;n<nOther;n++){
                                weight+= otherAssts[n]*otherMult[n];
                            }
                            otherWeights[m]=weight;
                            if(m!=(pow(2,nOther)-1)){
                                int max=0;
                                for(int n=0;n<nOther;n++){
                                    if(otherAssts[n]==0){
                                        max=n;
                                    }
                                }
                                otherAssts[max]=1;
                                for(int n=max+1;n<nOther;n++){
                                    otherAssts[n]=0;
                                }
                            }
                        }
                    }
                    // for each sampled draw, let us record the supported CPT(Y) given that Pa are the parents of Y
                    // cptsSupported will store the UNIQUE supported CPTs
                    // cptWeight will be their weights (# of 1s)
                    // cptsSupportCount will be the number of supporting draws 
                    std::vector<int> cptsSupported={};
                    std::vector<int> cptWeight={};
                    std::vector<int> cptsSupportCount={};
                    // for each sample:
                    for(int l=0;l<mcReps;l++){
                        // sample is the currently considered sample value set
                        std::vector<double> sample(nOut);
                        for(int m=0;m<nOut;m++){
                            sample[m]=FullSample[l*nOut + m];
                        }
                        // agreement gives, for each Pa assignment, which rule sample supports - 0 for 1>2, 1 for 2>1, 2 for both
                        std::vector<int> agreement(pow(2,nPa));
                        // paAsst cycles through the possible Pa assignments
                        std::vector<int> paAsst(nPa,0);
                        for(int m=0;m<pow(2,nPa);m++){
                            // calculate the fixed weight added to lexicographic posn number by the parental assignment:
                            int fixedWeight=0;
                            int counter=0;
                            for(int n=0;n<nvar;n++){
                                if(Pa[n]==1){
                                    fixedWeight+=paAsst[counter]*LexMult[n];
                                    counter+=1;
                                }
                            }
                            // group1 gives the lexicographic positions of the outcomes with this Pa asst and Y=0
                            // group2 gives the ones with Y=1
                            std::vector<int> group1(pow(2,nOther));
                            std::vector<int> group2(pow(2,nOther));
                            if(nOther==0){
                                group1[0]=fixedWeight;
                                group2[0]=fixedWeight+LexMult[col];
                            }
                            else{
                                for(int n=0;n<pow(2,nOther);n++){
                                    group1[n]=otherWeights[n]+fixedWeight;
                                    group2[n]=otherWeights[n]+fixedWeight+LexMult[col];
                                }
                            }
                            // value1/2 give the sum of the values of outcomes in group1/2 in the sample
                            double value1=0.0;
                            double value2=0.0;
                            for(int n=0;n<pow(2,nOther);n++){
                                value1+=sample[group1[n]];
                                value2+=sample[group2[n]];
                            }
                            //now that we know the sum values we can put which rule(s) the sample supports into the agreement vector
                            if(value1>value2){
                                agreement[m]=0;
                            }
                            else{
                                if(value2>value1){
                                    agreement[m]=1;
                                }
                                else{
                                    agreement[m]=2;
                                }
                            }
                            //go to the next Pa Asst
                            if(m<pow(2,nPa)-1){
                                int max=0;
                                for(int n=0;n<nPa;n++){
                                    if(paAsst[n]==0){
                                        max=n;
                                    }
                                }
                                paAsst[max]=1;
                                for(int n=max+1;n<nPa;n++){
                                    paAsst[n]=0;
                                }
                            }
                        }
                        //agreement now encodes which cpts are supported by sample
                        //let us now randomly select one by replacing each 2 entry by 0 or 1 randomly
                        struct timeval tv;
                        gettimeofday(&tv,0);
                        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                        std::default_random_engine eng(mySeed);
                        std::uniform_int_distribution<> UnifIntDist(0, 1);
                        for(int m=0;m<pow(2,nPa);m++){
                            if(agreement[m]==2){
                                agreement[m]=UnifIntDist(eng);
                            }
                        }
                        // let us add this to the list of supported CPTs if its new or add a support count if not
                        // first determine its weight
                        int weight=0;
                        for(int m=0;m<pow(2,nPa);m++){
                            weight+=agreement[m];
                        }
                        // how many CPTs are supported so far?
                        int nSupCpts=cptWeight.size();
                        if(nSupCpts==0){
                            // if no supported CPTs yet, add this one, its weight and 1 support count to the vectors respectively
                            for(int m=0;m<pow(2,nPa);m++){
                                cptsSupported.push_back(agreement[m]);
                            }
                            cptWeight.push_back(weight);
                            cptsSupportCount.push_back(1);
                        }
                        else{//If there are already some supported CPTs
                            //let us see if our CPT matches any of the others
                            bool NewCPT=true;
                            for(int m=0;m<nSupCpts;m++){
                                //for each existing supported CPT, see if it has same weight
                                if(cptWeight[m]==weight){
                                    //if they have equal weight, see if they're the same
                                    int ndiffs=0;
                                    for(int n=0;n<pow(2,nPa);n++){
                                        if(agreement[n]!=cptsSupported[(pow(2,nPa)*m)+n]){
                                           ndiffs+=1;
                                           break; 
                                        }
                                    }
                                    if(ndiffs==0){
                                        //if we have found it in the list then add a support count and stop looking
                                        cptsSupportCount[m]+=1;
                                        NewCPT=false;
                                        break;
                                    }
                                }
                            }
                            if(NewCPT){
                                //if we didnt find it in the list, add the CPT to the end of the list, add its weight, and assign a support count of 1
                                for(int m=0;m<pow(2,nPa);m++){
                                    cptsSupported.push_back(agreement[m]);
                                }
                                cptWeight.push_back(weight);
                                cptsSupportCount.push_back(1);
                            }
                        }
                    }
                    //cptsSupported now gives the unique CPTs supported by the Dirichlet sample
                    //cptsSupportCount gives their support counts, which is max?
                    int max=0;
                    int maxposn;
                    int nSupCpts=cptsSupportCount.size();
                    for(int l=0; l<nSupCpts; l++){
                        if(cptsSupportCount[l]>max){
                            max=cptsSupportCount[l];
                            maxposn=l;
                        }
                    }
                    //Let us identify ALL SupportedCPTs with this score and randomly select one to assign as the (potential) optimal CPT
                    int nBestCPTs=1;
                    vector<int> BestCPTs={maxposn};
                    for(int l=maxposn+1; l<nSupCpts; l++){
                        if(cptsSupportCount[l]==max){
                            BestCPTs.push_back(l);
                            nBestCPTs+=1;
                        }
                    }
                    struct timeval tv;
                    gettimeofday(&tv,0);
                    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                    std::default_random_engine eng(mySeed);
                    std::uniform_int_distribution<> UnifIntDist(0, nBestCPTs-1);
                    // ChosenCPT is the optimal CPT(Y) in the scenario where edge j->Y is removed (i.e. where Y has parents Pa)
                    int ChosenCPTPosn=BestCPTs[UnifIntDist(eng)];
                    // Delta value of j->Y is the table score of this optimal CPT/current Y score 
                    delta[j][col]=(double) cptsSupportCount[ChosenCPTPosn]/((double) mcReps * scores[col]);
                    // identify the optimal CPT and record in future CPTs
                    std::vector<int> cpt(pow(2,nPa+1));
                    int counter=0;
                    for(int l=0; l<pow(2,nPa); l++){
                        if(cptsSupported[(ChosenCPTPosn*pow(2,nPa))+l] == 0){
                            cpt[counter]=1;
                            cpt[counter+1]=2;
                        }
                        else{
                            cpt[counter]=2;
                            cpt[counter+1]=1;
                        }
                        counter+=2;
                    }
                    futureCPTs[j*nvar+col].resize(pow(2,nPa+1));
                    futureCPTs[j*nvar+col]=cpt;
                } 
                else{// if j->Y it is an edge to be added
                    // Let Pa be a 01 vector of the parents of Y and j
                    // nPa - number of parents in Pa
                    std::vector<int> Pa(nvar,0);
                    int nPa=0;
                    for(int m=0;m<nvar;m++){
                        if(m!=j){
                            if(cpnet.structure[m*nvar+col]==1){
                                Pa[m]=1;
                                nPa+=1;
                            }
                        }
                        else{
                            Pa[m]=1;
                            nPa+=1;
                        }
                    }
                    //We now calculate the score of Y in the scenario where it has parents Pa (i.e. when edge j->Y is added)
                    //This also gives us the optimal CPT in this scenario
                    //This calculation is done in the same way as all previous score calculations
                    //nother is the number of variables not Y or the variables in Pa (Ys parents plus j)
                    int nOther=nvar-1-nPa;
                    // otherWeights - list of all possible weights they can add to the lexicographic posn
                    std::vector<int> otherWeights(pow(2,nOther));
                    std::vector<int> otherAssts(nOther,0);
                    std::vector<int> otherMult(nOther);
                    if(nOther>0){
                        int counter=0;
                        for(int m=0;m<nvar;m++){
                            if(m!=col){
                                if(Pa[m]==0){
                                    otherMult[counter]=LexMult[m];
                                    counter+=1;
                                }
                            }
                        }
                        for(int m=0;m<pow(2,nOther);m++){
                            int weight=0;
                            for(int n=0;n<nOther;n++){
                                weight+= otherAssts[n]*otherMult[n];
                            }
                            otherWeights[m]=weight;
                            if(m!=(pow(2,nOther)-1)){
                                int max=0;
                                for(int n=0;n<nOther;n++){
                                    if(otherAssts[n]==0){
                                        max=n;
                                    }
                                }
                                otherAssts[max]=1;
                                for(int n=max+1;n<nOther;n++){
                                    otherAssts[n]=0;
                                }
                            }
                        }
                    }
                    // for each sampled draw, let us record the supported CPT(Y) in the scenario where Y has parents Pa
                    // cptsSupported will store the UNIQUE supported CPTs
                    // cptWeight will be their weights (# of 1s)
                    // cptsSupportCount will be the number of supporting draws 
                    std::vector<int> cptsSupported={};
                    std::vector<int> cptWeight={};
                    std::vector<int> cptsSupportCount={};
                    // for each sample:
                    for(int l=0;l<mcReps;l++){
                        // sample is the currently considered sample value set
                        std::vector<double> sample(nOut);
                        for(int m=0;m<nOut;m++){
                            sample[m]=FullSample[l*nOut + m];
                        }
                        // agreement gives, for each Pa asst, which rule sample supports - 0 for 1>2, 1 for 2>1, 2 for both
                        std::vector<int> agreement(pow(2,nPa));
                        // paAsst - cycles through the possible Pa assignments
                        std::vector<int> paAsst(nPa,0);
                        for(int m=0;m<pow(2,nPa);m++){
                            // calculate the fixed weight added to lexicographic posn number by the parental assignment:
                            int fixedWeight=0;
                            int counter=0;
                            for(int n=0;n<nvar;n++){
                                if(Pa[n]==1){
                                    fixedWeight+=paAsst[counter]*LexMult[n];
                                    counter+=1;
                                }
                            }
                            // group1 gives the positions of the outcomes with this parent asst and Y=0
                            // group2 gives the ones with Y=1
                            std::vector<int> group1(pow(2,nOther));
                            std::vector<int> group2(pow(2,nOther));
                            if(nOther==0){
                                group1[0]=fixedWeight;
                                group2[0]=fixedWeight+LexMult[col];
                            }
                            else{
                                for(int n=0;n<pow(2,nOther);n++){
                                    group1[n]=otherWeights[n]+fixedWeight;
                                    group2[n]=otherWeights[n]+fixedWeight+LexMult[col];
                                }
                            }
                            // value1/2 give the sum of the values of the outcomes in group1/2 in the sample
                            double value1=0.0;
                            double value2=0.0;
                            for(int n=0;n<pow(2,nOther);n++){
                                value1+=sample[group1[n]];
                                value2+=sample[group2[n]];
                            }
                            //now that we know the sum values we can put which rule(s) the sample supports into the agreement vector
                            if(value1>value2){
                                agreement[m]=0;
                            }
                            else{
                                if(value2>value1){
                                    agreement[m]=1;
                                }
                                else{
                                    agreement[m]=2;
                                }
                            }
                            //go to the next pa Asst
                            if(m<pow(2,nPa)-1){
                                int max=0;
                                for(int n=0;n<nPa;n++){
                                    if(paAsst[n]==0){
                                        max=n;
                                    }
                                }
                                paAsst[max]=1;
                                for(int n=max+1;n<nPa;n++){
                                    paAsst[n]=0;
                                }
                            }
                        }
                        //agreement now encodes which cpts are supported by sample
                        //let us now randomly select one by replacing each 2 entry by 0 or 1 randomly
                        struct timeval tv;
                        gettimeofday(&tv,0);
                        unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                        std::default_random_engine eng(mySeed);
                        std::uniform_int_distribution<> UnifIntDist(0, 1);
                        for(int m=0;m<pow(2,nPa);m++){
                            if(agreement[m]==2){
                                agreement[m]=UnifIntDist(eng);
                            }
                        }
                        // let us add this to the list of supported CPTs if its new or add a support count if not
                        // first determine its weight (# of 1s)
                        int weight=0;
                        for(int m=0;m<pow(2,nPa);m++){
                            weight+=agreement[m];
                        }
                        // how many CPTs are supported so far?
                        int nSupCpts=cptWeight.size();
                        if(nSupCpts==0){
                            // if no supported CPTs yet, add this one, its weight and 1 support count to the vectors respectively
                            for(int m=0;m<pow(2,nPa);m++){
                                cptsSupported.push_back(agreement[m]);
                            }
                            cptWeight.push_back(weight);
                            cptsSupportCount.push_back(1);
                        }
                        else{//If there are some supported CPTs already
                            //let us see if our CPT matches any of the others
                            bool NewCPT=true;
                            for(int m=0;m<nSupCpts;m++){
                                //for each existing supported CPT, see if it has same weight
                                if(cptWeight[m]==weight){
                                    //if they have equal weight, see if they're the same
                                    int ndiffs=0;
                                    for(int n=0;n<pow(2,nPa);n++){
                                        if(agreement[n]!=cptsSupported[(pow(2,nPa)*m)+n]){
                                           ndiffs+=1;
                                           break; 
                                        }
                                    }
                                    if(ndiffs==0){
                                        //if we have found our CPT in the list then add a support count and stop looking
                                        cptsSupportCount[m]+=1;
                                        NewCPT=false;
                                        break;
                                    }
                                }
                            }
                            if(NewCPT){
                                //if we didnt find the CPT in the list, add the new CPT to the end, add its weight, assign a support count of 1
                                for(int m=0;m<pow(2,nPa);m++){
                                    cptsSupported.push_back(agreement[m]);
                                }
                                cptWeight.push_back(weight);
                                cptsSupportCount.push_back(1);
                            }
                        }
                    }   
                    //cptsSupported now gives the unique CPTs supported by the Dirichlet sample
                    //cptsSupportCount gives the support counts, which is max?
                    int max=0;
                    int maxposn;
                    int nSupCpts=cptsSupportCount.size();
                    for(int l=0; l<nSupCpts; l++){
                        if(cptsSupportCount[l]>max){
                            max=cptsSupportCount[l];
                            maxposn=l;
                        }
                    }
                    //Let us identify ALL SupportedCPTs with maximal support and randomly select one to assign as the (potential) optimal CPT
                    int nBestCPTs=1;
                    vector<int> BestCPTs={maxposn};
                    for(int l=maxposn+1; l<nSupCpts; l++){
                        if(cptsSupportCount[l]==max){
                            BestCPTs.push_back(l);
                            nBestCPTs+=1;
                        }
                    }
                    struct timeval tv;
                    gettimeofday(&tv,0);
                    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
                    std::default_random_engine eng(mySeed);
                    std::uniform_int_distribution<> UnifIntDist(0, nBestCPTs-1);
                    // This optimal CPT is the maximally supported CPT(Y) in the scenario where Y has parents Pa (i.e. when j->Y has been added)
                    int ChosenCPTPosn=BestCPTs[UnifIntDist(eng)];
                    // Delta value of j->Y is table score of this optimal CPT/current Y score 
                    delta[j][col]=(double) cptsSupportCount[ChosenCPTPosn]/((double) mcReps * scores[col]);
                    // identify the optimal CPT and record this in future CPTs
                    // This is the CPT(Y) we will implement if j->Y is added
                    std::vector<int> cpt(pow(2,nPa+1));
                    int counter=0;
                    for(int l=0; l<pow(2,nPa); l++){
                        if(cptsSupported[(ChosenCPTPosn*pow(2,nPa))+l] == 0){
                            cpt[counter]=1;
                            cpt[counter+1]=2;
                        }
                        else{
                            cpt[counter]=2;
                            cpt[counter+1]=1;
                        }
                        counter+=2;
                    }
                    futureCPTs[j*nvar+col].resize(pow(2,nPa+1));
                    futureCPTs[j*nvar+col]=cpt;
                }
            }
        }
        

        //We have now updated the CP-net, the scores, the deltas, and the potentially optimal CPTs
        //We can now determine whether or not there are further valid changes we can make
        AdjMat=cpnet.structure;
        //Continue determines whether or not our learning process continues (i.e. if there are further valid changes)
        Continue=false;
        //ValidChanges gives the delta value of valid changes and -1 for invalid changes
        //For each edge j->k
        for(int j=0;j<nvar;j++){
            for(int k=0;k<nvar;k++){
            	// If j->k has -1 entry in cyclic, then changin this edge creates cycles so it is not valid
                if(cyclic[j][k]!=0){
                    ValidChanges[j*nvar +k]=-1;
                }
                else{
                	// If j->k can be added without creating cycles
                	// Determine whether it is an edge for removal or addition (i.e. if it is/is not already in structure) and then evaluate whether the delta value of j->k is above the appropriate threshold
                	// If it is above, then it is a valid edge change (add delta value to ValidChanges). If not, it is not a valid change, add -1 entry
                    if(AdjMat[j*nvar +k]==0){
                        if(delta[j][k]>(1+signif/100)){
                            ValidChanges[j*nvar +k]=delta[j][k];
                            Continue=true;
                        }
                        else{
                            ValidChanges[j*nvar +k]=-1;
                        }
                    }
                    else{
                        if(delta[j][k]>(1/(1+signif/100))){
                            ValidChanges[j*nvar +k]=delta[j][k];
                            Continue=true;
                        }
                        else{
                            ValidChanges[j*nvar +k]=-1;
                        }
                    }
                }
            }
        }

        //Continue is now true if and only if there is at least one valid edge change.
        //If there are more edge changes we can make, we repeat this while loop, implementing the best possible edge change
        //If there are no more edge changes possible our learning process is over and we exit the loop

        vector<int>breaks=cpnet.cptBreaks;
        vector<int> entries=cpnet.cptEntries;
    }
    //Our learning procedure  is now over
    //Add the final, learned CP-net to our learning sequence (essential as if no steps were taken it is no yet recorded)
    LearnedCpns.push_back (cpnet);

    //Store learning sequence and final scores in output
    struct LearningOutput output;
    output.LearnedCpnSeq=LearnedCpns;
    output.FinalScores=scores;
    //Evaluate the time elapsed and store this in output
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double TimeJump=elapsed.count();
    output.TimeElapsed=TimeJump;
    
    //Write out the sequence of learned CP-nets to the specified file
    int nCpns=LearnedCpns.size();
    for(int i=0;i<nCpns;i++){
        writeCpn(CpnSeqFilename, LearnedCpns[i]);
    }
    // Write the final scores to this file
    ofstream myfile (CpnSeqFilename, ios::app);
    for(int i=0;i<nvar;i++){
        myfile << scores[i];
        if(i<nvar-1){
            myfile << ",";
        }
    }
    myfile << "\n";
    // Write out the time jump
    myfile << TimeJump;
    myfile << "\n";
    myfile.close();

    //We write out the results so they can be saved and further processed via testing etc
    //We also return the restults so that the learning result can be used directly within other functions.
    return output;
}
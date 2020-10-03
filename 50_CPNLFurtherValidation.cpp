#include "50ValidationMeasures.h"

vector<double> EntailmentAgreement(cpn NLearned, vector<vector<int>> DQ1, vector<vector<int>> DQ2){
    //This function will take the learned CP-net and a list of preferences entailed by the original/true CP-net
    // Entailed preferences are a list of N pairs o1>o2
    // DQ1 is a vector of N outcomes, the o1 outcomes
    // DQ2 contains the o2 outcomes
    // NLearned is the learned CP-net
    //We return 3 values (regarding the list of entailed preferences)
    //-the number that the learned CP-net entails
    //-the number that the learned CP-net entails the opposite direction
    //-the number of cases where the learned CP-net entails neither preference

    // Extract number of variables and domain sizes from NLearned
    vector<int> Doms=NLearned.domains;
    int nvar=Doms.size();
    //start by setting the three values to 0
    double Agree=0.0;
    double Disagree=0.0;
    double Indecisive=0.0;
    //nPrefs - number of input entailed preferences
    int nPrefs=DQ1.size();
    //add 1 to CPN breaks to make the CP-nets compatible with dominance testing function (left over issue from R to C++ translation)
    vector<int> breaks=NLearned.cptBreaks;
    for(int i=0;i<breaks.size();i++){
      breaks[i]+=1;
    }
    //For each preference, determine the direction entailed by NLearn
    for(int i=0;i<nPrefs;i++){
        //the preference we are considering:
        vector<int> outcome1=DQ1[i];
        vector<int> outcome2=DQ2[i];
        //does CPNL entail the correct direction? (Use dominance testing function to answer)
        List DTResult=RSDQRankPriority(NLearned.structure, NLearned.domains, NLearned.cptEntries, breaks, outcome1, outcome2);
        //if so, add a count to agree
        if(DTResult.result){
            Agree+=1;
        }
        else{
            //does CPNL entail the wrong direction?
            DTResult=RSDQRankPriority(NLearned.structure, NLearned.domains, NLearned.cptEntries, breaks, outcome2, outcome1);
            //if so, add a count to disagree
            if(DTResult.result){
                Disagree+=1;
            }
            else{
                //if neither direction is entailed, add a count to indecisive
                Indecisive+=1;
            }
        }
    }
    //Return these counts
    vector<double> Counts(3);
    Counts[0]=Agree;
    Counts[1]=Disagree;
    Counts[2]=Indecisive;
    return Counts;
}

int PrefGraphMisorientation(cpn NOriginal, cpn NLearned){
    //This function will take 2 binary CPNs of the same size and count the number of edges in their preference graphs that
    //have different orientations
    int WrongEdges=0;
    //obtain the structures (Ag - original, Al - learned)
    vector<int> Ag=NOriginal.structure;
    vector<int> Al=NLearned.structure;
    //obtain the number of variables
    vector<int> Doms=NOriginal.domains;
    int nvar=Doms.size();
    //obtain the breaks and entries vectors of both CP-nets
    vector<int> Breaksg=NOriginal.cptBreaks;
    vector<int> Breaksl=NLearned.cptBreaks;
    vector<int> Entriesg=NOriginal.cptEntries;
    vector<int> Entriesl=NLearned.cptEntries;
    //for each variable identify the number flips of that variable which are oriented differently in the preference graphs of NOriginal and NLearned
    for(int i=0; i<nvar; i++){
        //obtain parent sets of variable i in both CP-nets
        vector<int> Pag(nvar,0);
        vector<int> Pal(nvar,0);
        //Intersect - variables that are a parent of variable i in both CP-nets
        vector<int> Intersect(nvar,0);
        int nPag=0;
        int nPal=0;
        int nInt=0;
        for(int j=0;j<nvar;j++){
            int both=0;
            if(Ag[j*nvar +i]==1){
                Pag[j]=1;
                nPag+=1;
                both+=1;
            }
            if(Al[j*nvar +i]==1){
                Pal[j]=1;
                nPal+=1;
                both+=1;
            }
            if(both==2){
                Intersect[j]=1;
                nInt+=1;
            }
        }
        //Let us find the CPTs for variable j (which we will call X) in each CP-net
        vector<int> CPTXg(Breaksg[i+1]-Breaksg[i]);
        //CPTleng - length of CPT(X) in NOriginal 
        int CPTleng=Breaksg[i+1]-Breaksg[i];
        //Extract CPT(X) from NOriginal entries vector
        for(int j=0;j<Breaksg[i+1]-Breaksg[i];j++){
            CPTXg[j]=Entriesg[Breaksg[i]+j];
        }
        //CPTlenl - length of CPT(X) in NLearned 
        int CPTlenl=Breaksl[i+1]-Breaksl[i];
        vector<int> CPTXl(Breaksl[i+1]-Breaksl[i]);
        //Extract CPT(X) from NLearned entries vector
        for(int j=0;j<Breaksl[i+1]-Breaksl[i];j++){
            CPTXl[j]=Entriesl[Breaksl[i]+j];
        }
        if(nInt>0){
            //if there IS some intersection in parents of X
            if((nInt==nPag)&&(nInt==nPal)){
                //if the parent set of X is the same in both:
                //Number of variables other than variable and parents
                int nOther=nvar-1-nPag;
                //The mis-oriented X flips will be 2^nOther * #CPT(X) rules where NOriginal and NLearned differ
                //MisOrientedRules- #CPT(X) rules where NOriginal and NLearned differ
                int MisOrientedRules=0;
                for(int j=0;j<(CPTleng/2);j++){
                    if(CPTXg[2*j]!=CPTXl[2*j]){
                        MisOrientedRules+=1;
                    }
                }
                WrongEdges+=MisOrientedRules*pow(2,nOther);
            }
            else{
                //if parent sets of X are not equal OR disjoint:
                //number of variables not the variable or either parent set
                int nOther=nvar-1-nPag-nPal+nInt;
                //Lexicographic multipliers for a parent assignment are the numbers you multiply the assignment by (then sum) to find the lexicographic order number of a given parent assignment
                //IntLexMultg - The lexicographic multipliers for the intersection variables, considered as part of a parent assignment (of Pa(X)) in NOriginal
                //OtherLexMultg - The lexicographic multipliers for the parents of X not in intersection in NOriginal
                // Similarly for IntLexMultl and OtherLexMultl
                vector<int> IntLexMultg(nInt);
                vector<int> IntLexMultl(nInt);
                vector<int> OtherLexMultg(nPag-nInt);
                vector<int> OtherLexMultl(nPal-nInt);
                int posncounter1=nInt-1;
                int posncounter2=nPag-nInt-1;
                int valuecounter=0;
                for(int j=nvar-1;j>=0;j--){
                    if(Pag[j]==1){
                        if(Intersect[j]==1){
                            IntLexMultg[posncounter1]=valuecounter;
                            posncounter1-=1;
                            valuecounter+=1;
                        }
                        else{
                            OtherLexMultg[posncounter2]=valuecounter;
                            posncounter2-=1;
                            valuecounter+=1;
                        }
                    }
                }
                posncounter1=nInt-1;
                posncounter2=nPal-nInt-1;
                valuecounter=0;
                for(int j=nvar-1;j>=0;j--){
                    if(Pal[j]==1){
                        if(Intersect[j]==1){
                            IntLexMultl[posncounter1]=valuecounter;
                            posncounter1-=1;
                            valuecounter+=1;
                        }
                        else{
                            OtherLexMultl[posncounter2]=valuecounter;
                            posncounter2-=1;
                            valuecounter+=1;
                        }
                    }
                }
                for(int j=0; j<nInt;j++){
                    IntLexMultg[j]=pow(2,IntLexMultg[j]);
                    IntLexMultl[j]=pow(2,IntLexMultl[j]);
                }
                for(int j=0; j<nPag-nInt;j++){
                    OtherLexMultg[j]=pow(2,OtherLexMultg[j]);
                }
                for(int j=0; j<nPal-nInt;j++){
                    OtherLexMultl[j]=pow(2,OtherLexMultl[j]);
                }
                //for every assignment of values to the intersection:
                vector<int> IntAsst(nInt,0);
                for(int j=0;j<pow(2,nInt);j++){
                    //calculate the fixed weight contributed to parental assignment lexicographic position (for both NOriginal and NLearned) by this intersection assignment
                    int FixedWeightg=0;
                    int FixedWeightl=0;
                    for(int k=0;k<nInt;k++){
                        if(IntAsst[k]==1){
                            FixedWeightg+=IntLexMultg[k];
                            FixedWeightl+=IntLexMultl[k];
                        }
                    }
                    //count the number of rules in each direction with this asst to intersection
                    //nDirection1g - number of 1>2 rules in CPT(X) in NOriginal where parent assinment sets the intersection to IntAsst
                    //nDirection2g - number of 2>1 rules in CPT(X) in NOriginal where parent assinment sets the intersection to IntAsst
                    //nDirection1l and nDirection2l defined similarly for NLearned
                    int nDirection1g=0;
                    int nDirection1l=0;
                    int nDirection2g=0;
                    int nDirection2l=0;
                    //for each possible asssignment to the non-intersection parents in NOriginal, determine the CPT rule and add a count to the appropriate variable
                    if(nPag-nInt==0){
                        //if Pa(X)=intersection in NOriginal, there is only one rule to evaluate (i.e that corresponding to parents=IntAsst)
                        if(CPTXg[2*FixedWeightg]==1){
                            nDirection1g+=1;
                        }
                        else{
                            nDirection2g+=1;
                        }
                    }
                    else{
                        //OtherAsst cycles through all possible assignments. This, together with IntAsst fully specifies the pa(X) assignment
                        //For each OtherAsst we determine the preference rule and add a count to the appropriate direction
                        vector<int> OtherAsst(nPag-nInt,0);
                        for(int k=0;k<pow(2,nPag-nInt);k++){
                            //parentposn is the lexicographic position of this parental assignment
                            int parentposn=FixedWeightg;
                            for(int l=0;l<nPag-nInt;l++){
                                parentposn+=OtherAsst[l]*OtherLexMultg[l];
                            }
                            //if the rule under this parent assignment is 1>2, we add a count to direction 1
                            if(CPTXg[2*parentposn]==1){
                                nDirection1g+=1;
                            }
                            else{//if the rule under this parent assignment is 2>1, we add a count to direction 2
                                nDirection2g+=1;
                            }
                            //Move to the next assignment of OtherAsst
                            if(k<pow(2,nPag-nInt)-1){
                                int max=0;
                                for(int l=0;l<nPag-nInt;l++){
                                    if(OtherAsst[l]==0){
                                        max=l;
                                    }
                                }
                                OtherAsst[max]=1;
                                for(int l=max+1;l<nPag-nInt;l++){
                                    OtherAsst[l]=0;
                                }
                            }
                        }
                    }
                    //We now add counts to nDirection1l and nDirection2l in the same way, now using NLearned
                    if(nPal-nInt==0){
                        if(CPTXl[2*FixedWeightl]==1){
                            nDirection1l+=1;
                        }
                        else{
                            nDirection2l+=1;
                        }
                    }
                    else{
                        vector<int> OtherAsst(nPal-nInt,0);
                        for(int k=0;k<pow(2,nPal-nInt);k++){
                            int parentposn=FixedWeightl;
                            for(int l=0;l<nPal-nInt;l++){
                                parentposn+=OtherAsst[l]*OtherLexMultl[l];
                            }
                            if(CPTXl[2*parentposn]==1){
                                nDirection1l+=1;
                            }
                            else{
                                nDirection2l+=1;
                            }
                            if(k<pow(2,nPal-nInt)-1){
                                int max=0;
                                for(int l=0;l<nPal-nInt;l++){
                                    if(OtherAsst[l]==0){
                                        max=l;
                                    }
                                }
                                OtherAsst[max]=1;
                                for(int l=max+1;l<nPal-nInt;l++){
                                    OtherAsst[l]=0;
                                }
                            }
                        }
                    }
                    //Now that we have these values, we can calculate the number of X-flips where Intersection=IntAsst and the two preference graphs have different orientation
                    WrongEdges+=(nDirection1g*nDirection2l+nDirection2g*nDirection1l)*pow(2,nOther);
                    //Now move to the next IntAsst
                    if(j<pow(2,nInt)-1){
                        int max=0;
                        for(int k=0;k<nInt;k++){
                            if(IntAsst[k]==0){
                                max=k;
                            }
                        }
                        IntAsst[max]=1;
                        for(int k=max+1;k<nInt;k++){
                            IntAsst[k]=0;
                        }
                    }
                }
                //WrongEdges has now added up the number of misoriented (not matching) X flips for every possible intersection assignment, meaning that WrongEdges is now the total number of misoriented X flips
            }
        }
        else{
            //if there is NO intersection between parent sets:
            //Number of variables not X or parents
            int nOther=nvar-1-nPag-nPal;
            //count the number of rules in each direction for each CPT
            //nDirection1g - number of rules 1>2 in CPT(X) for NOriginal. The other terms are similarly defined.
            int nDirection1g=0;
            int nDirection1l=0;
            int nDirection2g=0;
            int nDirection2l=0;
            for(int j=0;j<CPTleng/2;j++){
                if(CPTXg[2*j]==1){
                    nDirection1g+=1;
                }
                else{
                    nDirection2g+=1;
                }
            }
            for(int j=0;j<CPTlenl/2;j++){
                if(CPTXl[2*j]==1){
                    nDirection1l+=1;
                }
                else{
                    nDirection2l+=1;
                }
            }
            //given these values, we can calculate the number of X flips oriended differently in NOriginal and NLearned as follows:
            WrongEdges+=(nDirection1g*nDirection2l+nDirection2g*nDirection1l)*pow(2,nOther);
        }
    }
    //WrongEdges now has summed the number of misoriented (not matching) flips of every variable
    //Thus, WrongEdges is now the total number of misoriented edges in the preference graphs 
    return WrongEdges;
}

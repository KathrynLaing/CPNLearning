#include "domtest.h"
#include "cpnInput.h"


/* Take CP- net, variable X and parent assignment as inputs, output preference row corresponding to parental assignment */
/* No Parents - parent asst input should just be the number 0*/
/*variables and variable assignments (including CPT entries) indexed from 0*/
vector<int> CPTRow(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, int x, vector<int> PA){
  //extract CPT(X) from entries
  int x0=BREAK[(x-1)];
  int x1=BREAK[(x)];
  vector<int> CPT((x1-x0));
  for(int i=0; i<CPT.size();i++){
            CPT[i]=ENTRIES[(x0+i-1)];
  }
  /*CPT is the CPT of variable x - the entries only, flattened out into a vector*/
  int nPa=PA.size();
  /*If there are no parents of X*/
  if(PA[0]==0){
    //CPT only has one parent assignment and only one row  so return this
    return CPT;
  }
  else{
    //PaDom is the domain sizes of the parents of X in order
    vector<int> PaDom(nPa);
    int j=0;
    //Note N.Size is the number of variables
    for(int i=0;i<N.size();i++){
      //If i is a parent of x:
      if(A[i*(N.size())+(x-1)]!=0){
        PaDom[j]=N[i];
        j+=1;
      }
    }
    //nPaAsst is the number associated with the input parent assignment when parent asst are enumerated lexicographically
    int nPaAsst;
    if(PA.size()==1){
      //If there is 1 parent only
      nPaAsst=PA[0];
    }
    else{//if there are multiple parents:
      nPaAsst=1;
      //iter is the lexicographic multiplier vector that we use to change between parental assignments and their lexicographic enumeration
      vector<int> iter(nPa);
      for(int i=0;i<nPa;i++){
        int prod=1;
        for(int j=i+1;j<nPa;j++){
          prod=prod*PaDom[j];
        }
        iter[i]=prod;
      }
      for(int i=0;i<nPa;i++){
        if(PA[i]!=1){
          int a=PA[i]-1;
          nPaAsst+=(a*iter[i]);
        }
      }
    }
    nPaAsst-=1;
    //dom - domain size of X
    int dom=N[(x-1)];
    //Extract the nPaAsst-th row of CPT(X). This will be the preference row corresponding to the input parent assignment
    vector<int> cptRow(dom);
    for(int i=0;i<dom;i++){
      cptRow[i]=CPT[(nPaAsst*dom)+i];
    }
    return cptRow;
  }
  /*We have returned the preference row of CPT(X) corresp to given parental asst */
}



double Rank(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, vector<int> o){
  // Input a CPN and an outcome o
  //Output r(o)
  //Start by setting rank to 0
  double rank=0.0;
  //AF - vector of ancestral factors
  //Note N.Size is the number of variables
  std::vector<double> AF(N.size());
  //AncSizes - vector of #ancestors for each variable
  vector<int> AncSizes(N.size());
  //For each variable: 
  for(int i=0;i<N.size();i++){
    //calculate the ancestor setof variable i
    //Start by adding parents and then repeatedly add the parents of every variable in the set until there is no change
    vector<int> Anc(N.size(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<N.size();j++){
      if(A[j*(N.size())+i]==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<N.size();j++){
        if(Anc[j]==1){
          for(int k=0;k<N.size();k++){
            if(A[k*(N.size())+j]==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<N.size();j++){
        AncSum+=Anc[j];
      }
    }
    //Anc is now a 0/1 vector giving the ancestors of variable i
    //AncSum - number of ancestors (non-zero entries in Anc)
    AncSizes[i]=AncSum;
    //Calculate the ancestral factor of variable i (AFi) by dividing 1 by the domain size of every ancestor of i
    double AFi=1.0;
    for(int j=0;j<N.size();j++){
      if(Anc[j]==1){
        AFi=AFi/N[j];
      }
    }
    AF[i]=AFi;
  }
  //DP- vector giving the number of descendent paths for each variable
  vector<int> DP(N.size(),-1);
  int counter=0;
  int index;
  //until all DP values are calculated:
  while(counter<N.size()){
    //identify a variable (index) with maximal ancestor set which has not yet had descendent paths calculated
    //by definition all of the children of index have already had DP value calculated.
    index=-1;
    int ancsize=-1;
    for(int i=0;i<N.size();i++){
      if(DP[i]==(-1)){
        if(AncSizes[i]>ancsize){
          index=i;
          ancsize=AncSizes[i];
        }
      }
    }
    //calculate DP of index by summing (DP+1) over all children of index
    int dp=0;
    for(int i=0;i<N.size();i++){
      if(A[index*(N.size())+i]==1){
        dp+=(DP[i]+1);
      }
    }
    DP[index]=dp;
    //increase counter of #DP values calculated
    counter+=1;
  }
  //For each variable, calculate the weight contributed to rank
  for(int i=0;i<N.size();i++){
    //First, we must determine the preference position of the value taken by variable i in o
    //Parents 0/1 vector of the parents of i
    //counter - number of parents of i
    vector<int> Parents(N.size(),0);
    counter=0;
    for(int j=0;j<N.size();j++){
      if(A[j*(N.size())+i]==1){
        Parents[j]=1;
        counter+=1;
      }
    }
    vector<int> pref;
    if(counter==0){
      //if theres no parents:
      vector<int> null(1,0);
      //pref - the CPT of i (only one row as no parents)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), null);
    }
    else{
      // if i does have parents
      //PA values assigned to parents of i in o
      vector<int> PA(counter);
      counter=0;
      for(int j=0;j<N.size();j++){
        if(Parents[j]==1){
          PA[counter]=o[j];
          counter+=1;
        }
      }
      //pref - The row of CPT(i) corresponding to the parent assigment given in o (to the parents of i)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PA);
    }
    //pp - position of preference term in rank formula
    //This is calculable with |Dom(i)| because we can extract the preference position of i in o from pref (as this is the preference rule over i corresponding to the parental assignment in o)
    double pp=((N[i]-pref[o[i]-1]+1));
    pp=pp/(N[i]);
    //we can now calculate the weight contributed to the rank by variable i as follows (using the rank formula)
    rank+=AF[i]*(DP[i]+1)*pp;
  }
  //rank is now r(o) as we started at 0 and summed the weights contributed by all variables.
  return rank;
}

std::vector<double> RankDifferences(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK){
  //Input - CP-net
  //Output - vector of the least rank improvements
  //We will store the least rank improvement in Diff
  std::vector<double> Diff(N.size());
  //AF - will be a vector of ancestral factors
  std::vector<double> AF(N.size());
  //AncSizes - will store the number of ancestors of each variable
  vector<int> AncSizes(N.size());
  //For each variable, calculate the ancestral factor:
  for(int i=0;i<N.size();i++){
    //First calculate the ancestor set (Anc) of variable i as in the rank function:
    vector<int> Anc(N.size(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<N.size();j++){
      if(A[j*(N.size())+i]==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<N.size();j++){
        if(Anc[j]==1){
          for(int k=0;k<N.size();k++){
            if(A[k*(N.size())+j]==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<N.size();j++){
        AncSum+=Anc[j];
      }
    }
    //AncSum is now the #ancestors of variable i
    AncSizes[i]=AncSum;
    //Calculate ancestral factor of i (AFi) as we did in the rank function:
    double AFi=1.0;
    for(int j=0;j<N.size();j++){
      if(Anc[j]==1){
        AFi=AFi/N[j];
      }
    }
    AF[i]=AFi;
  }
  //AF now contains ancestral factors and AncSizes contains #ancestors for each variable
  //DP - number of descendent paths for each variable, calculated in the same way as in the rank function
  vector<int> DP(N.size(),-1);
  int counter=0;
  int index;
  while(counter<N.size()){
    index=-1;
    int ancsize=-1;
    for(int i=0;i<N.size();i++){
      if(DP[i]==(-1)){
        if(AncSizes[i]>ancsize){
          index=i;
          ancsize=AncSizes[i];
        }
      }
    }
    int dp=0;
    for(int i=0;i<N.size();i++){
      if(A[index*(N.size())+i]==1){
        dp+=(DP[i]+1);
      }
    }
    DP[index]=dp;
    counter+=1;
  }
  //For each variable, calculate the least rank improvement term
  for(int i=0;i<N.size();i++){
    //difference - least rank improvement of variable i
    double difference=AF[i]*(DP[i]+1);
    difference=difference/(N[i]);
    double ChildTerms=0.0;
    //for each child of i, calculate the deduction term
    for(int j=0;j<N.size();j++){
      if(A[i*(N.size())+j]==1){
        double term=AF[j]*(DP[j]+1)*(N[j]-1);
        term=term/(N[j]);
        ChildTerms+=term;
      }
    }
    difference-=ChildTerms;
    //difference is now the least rank improvement of variable i, calculated according to the given formula
    Diff[i]=difference;
  }
  //Diff - vector of least rank improvement terms for each variable.
  return Diff;
}


List RSDQRankPriority(vector<int> A, vector<int> N, vector<int> ENTRIES, vector<int> BREAK, vector<int> o1, vector<int> o2) {
  //Function takes a CP-net and a dominance query (is o1>o2?) and answers it using rank + diff pruning with rank prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note N.size is the number of variables
  for(int i=0;i<N.size();i++){
    Omega*=N[i];
  }
  int n=N.size();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R;
    R.result=false;
    R.OT=0;
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R;
    R.result=false;
    R.OT=0;
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(N.size());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //AncSizes - vector of #ancestors for each variable
  vector<int> AncSizes(N.size());
  for(int i=0;i<N.size();i++){
    //For variable i, calculate the ancestor set in the same way as rank function
    vector<int> Anc(N.size(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<N.size();j++){
      if(A[j*(N.size())+i]==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<N.size();j++){
        if(Anc[j]==1){
          for(int k=0;k<N.size();k++){
            if(A[k*(N.size())+j]==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<N.size();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(N.size());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<N.size();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives rank value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=rank2;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcome (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3
    vector<int> o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        vector<int> Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A[j*(N.size())+i]==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parent assignment in o3
        vector<int> pref;
        vector<int> PaAsst(NPa);
        if(NPa==0){
          vector<int> null(1,0);
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), null);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          vector<int> flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            vector<int> o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
                //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
                //return true and the number of outcomes in the search tree
                long long int Count=SearchTree.size();
                List R;
                R.result=true;
                R.OT=Count;
                return R;
            }
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree:
                //add o4 to the search Tree and put its rank into SearchTreeRanks
                SearchTree.push_back (lex4);
                SearchTreeRanks.push_back (rank4);
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R;
  R.result=false;
  R.OT=Count;
  return R;
}


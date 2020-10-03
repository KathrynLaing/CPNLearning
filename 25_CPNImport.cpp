#include "cpnInput.h"

struct cpn readCpn(string filename, int n){
    //Takes a filepath and number of variables and reads in to CP-net structure
    //inputN will be the CP-net we build
    struct cpn inputN;
    // n is the number of variables
    inputN.size=n;
    // read in the adjacency matrix of CP-net structure row by row
    // Enter this into vector A where each row is listed end to end
    std::vector<int> A(n*n);
    int counter=0;
    string inputnumber;
    char cNum[3*n];
    ifstream infile;
    infile.open (filename, ifstream::in);
    for(int i=0;i<n;i++){
        infile.getline(cNum, 256);
        stringstream input_string(cNum);
        for(int j=0;j<n;j++){
            getline(input_string,inputnumber,',');
            A[counter]=atoi(inputnumber.c_str());
            counter++;
        }
    }
    inputN.structure=A;
    //read the domain sizes from the next line and pass to vector N
    std::vector<int> N(n);
    infile.getline(cNum, 256);
    stringstream input_string(cNum);
    for(int i=0; i<n; i++){
        getline(input_string,inputnumber,',');
        N[i]=atoi(inputnumber.c_str());
    }
    inputN.domains=N;
    //Use structure to calculate the CPT sizes, simultaneously creating cpt Breaks vector
    //NOTE OUR BREAKS WILL START AT 0  DUE TO C++ INDEXING WHEREAS THE BREAKS IN THE FILE WILL START AT 1
    std::vector<int> Breaks(n+1);
    //First entry of breaks is 0 - the starting point of CPT1
    Breaks[0]=0;
    int cptSize=0;
    //for each variable X
    for(int i=0; i<n; i++){
        //CPT(X) size starts at |Dom(X)| and is then multiplied by the domain size of each parent. 
        int cpti=N[i];
        for(int j=0; j<n; j++){
            if(A[j*n+i]==1){
            cpti*=N[j];    
            }
        }
        cptSize+=cpti;
        //cpti is the length of  CPT(X) in the entries vector. Use this to calculate where the next CPT starts and input to breaks
        Breaks[i+1]=Breaks[i]+cpti;
    }
    //cptSize is now the sum of all the CPT lengths and, thus, is the lenght of the entries vector
    //Breaks is now fully determined so we can add this to inputN
    inputN.cptBreaks=Breaks;
    //read in entries from next line in file and pass to Entries vector
    std::vector<int> Entries(cptSize);
    char line[3*cptSize];
    infile.getline(line, 3*cptSize);
    stringstream input_string2(line);
    for(int i=0; i<cptSize; i++){
        getline(input_string2,inputnumber,',');
        Entries[i]=atoi(inputnumber.c_str());
    }
    inputN.cptEntries=Entries;
    infile.close();
    //inputN is now fully specified from the data in the file
    return inputN;
}

vector<cpn> readMultCpn(string filename, int n){
    //This function will take a text file that may have one or more CPNs (one after the other) and convert them into a vector of CPNs
    //All CP-nets in file must have n variables
    //First, we find the number of lines in the file
    int Lcount=0;
    string line;
    ifstream LineFile(filename);
    while(getline(LineFile,line)){
        Lcount++;
    }
    //nCpns is the number of cpns in the file (each CP-net takes up n+3 lines of text)
    int nCpns=Lcount/(n+3);
    //cpnSeq will be the CPN vector we output
    vector<cpn> CpnSeq(nCpns);
    //open the file again
    ifstream infile(filename);
    //For each CP-net in the document we read it in in the same way as the previous funcion and then store the CP-net in CpnSeq
    for(int i=0;i<nCpns;i++){
        //create the next CPN
        struct cpn inputN;
        // n is the number of variables
        inputN.size=n;
        // read in CP-net structure and record in vector A
        std::vector<int> A(n*n);
        int counter=0;
        string inputnumber;
        char dummy;
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                infile >> A[counter];
                if(j<n-1){
                    infile >> dummy;
                }
                counter++;
            }
        }
        inputN.structure=A;
        //read the domain sizes from the next line and pass to vector N
        std::vector<int> N(n);
        for(int i=0; i<n; i++){
            infile >> N[i];
            if(i<n-1){
                infile >> dummy;
            }
        }
        inputN.domains=N;
        //Use structure to calculate the CPT sizes, simultaneously creating cpt Breaks vector, as we did above
        //NOTE OUR BREAKS AGAIN WILL START AT 0  DUE TO C++ INDEXING 
        std::vector<int> Breaks(n+1);
        Breaks[0]=0;
        int cptSize=0;
        for(int i=0; i<n; i++){
            int cpti=N[i];
            for(int j=0; j<n; j++){
                if(A[j*n+i]==1){
                cpti*=N[j];    
                }
            }
            cptSize+=cpti;
            Breaks[i+1]=Breaks[i]+cpti;
        }
        //Breaks is now fully determined so we can add this to inputN
        inputN.cptBreaks=Breaks;
        //read in entries from next line in file and pass to Entries vector
        std::vector<int> Entries(cptSize);
        for(int i=0; i<cptSize; i++){
            infile >> Entries[i];
            if(i<cptSize-1){
                infile >> dummy;
            }
        }
        inputN.cptEntries=Entries;
        //Read the breaks to a dummy variable
        //This essentially "skips" the line (as we already have the breaks vector) so that we are in the right place in the file to start reading the next vector
        int x;
        for(int i=0;i<n+1;i++){
            infile >> x;
            if(i<n){
                infile >> dummy;
            }
        }
        //Store the created CP-net in CpnSeq
        CpnSeq[i]=inputN;
    }
    //All CP-nets in the file have now been read in and stored in CpnSeq
    //Close the file
    infile.close();
    return CpnSeq;
}

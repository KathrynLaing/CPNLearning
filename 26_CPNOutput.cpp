#include "cpnOutput.h"
#include "cpnInput.h"

int writeCpn(string filename, cpn N){
    //This function takes a filename and a CP-net structure and writes it to the file in such a way it can be read in by the CPN Input function
    //If the file does not exist already it creates
    //If it does, it appends the CPN to the end
    //First, extract the number of variables and the structure from N
    int nvar=N.size;
    vector<int> AdjMat=N.structure;
    ofstream myfile (filename, ios::app);
    //write the srtucture to file, formatted as an nxn adjacency matrix
    for(int i=0;i<nvar;i++){
        for(int j=0;j<nvar;j++){
            myfile << AdjMat[i*nvar + j];
            if(j!=(nvar-1)){
                myfile << ",";
            }
        }
        myfile << "\n";
    }
    //write out the domain vector (n length vector giving variable domain sizes) on a new line
    vector<int> doms=N.domains;
    for(int i=0;i<nvar;i++){
        myfile << doms[i];
            if(i!=(nvar-1)){
                myfile << ",";
            }
    }
    myfile << "\n";
    //write out the CPT entries vector on a new line
    vector<int> entries=N.cptEntries;
    int cptSize=entries.size();
    for(int i=0;i<cptSize;i++){
        myfile << entries[i];
            if(i!=(cptSize-1)){
                myfile << ",";
            }
    }
    myfile << "\n";
    //write out the CPT breaks  vector on a new line
    //Note, we add 1 to the breaks so that they're indexed from 1 (not 0 as we use here). This is so output format matches up with the R output format
    //This enables us to use the same read-in function for all CP-nets.
    vector<int> breaks=N.cptBreaks;
    for(int i=0;i<nvar+1;i++){
        myfile << breaks[i]+1;
            if(i!=(nvar)){
                myfile << ",";
            }
    }
    myfile << "\n";
    myfile.close();
    return 0;
}

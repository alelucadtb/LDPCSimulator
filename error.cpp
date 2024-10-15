#include "error.h"

Error::Error(){

}

Error::~Error(){
    
}  

double Error::calculateError(std::vector<int>& received, Word& original){
    double error = 0.0;
    for(int i = 0; i < received.size(); i++){
        if(received[i] != original.get(i)){
            error++;
        }
    }
    return error/(double)received.size();
}

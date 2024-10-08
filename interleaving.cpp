#include "interleaving.h"

Interleaving::Interleaving(){

}

Interleaving::~Interleaving() = default;

std::vector<std::vector<double>> Interleaving::interleaving(std::vector<std::vector<int>> matrix){
    std::vector<std::vector<double>> result;
    std::vector<double> resultTmp;
    // For a fixed column of the matrix
    for(int i = 0; i < matrix[0].size(); i++){
        // For all the rows of the matrix
        for(int j = 0; j < matrix.size(); j++){
            resultTmp.push_back((double)matrix[j][i]);
        }
        result.push_back(resultTmp);
        resultTmp.clear();
    }
    return result;
}

std::vector<std::vector<double>> Interleaving::deinterleaving(std::vector<std::vector<double>> matrix){
    std::vector<std::vector<double>> result;
    std::vector<double> resultTmp;
    for(int i = 0; i < matrix[0].size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            resultTmp.push_back(matrix[j][i]); 
        }
        result.push_back(resultTmp);
        resultTmp.clear();
    }
    return result;
}
// C++ Program to Implement Adjacency List
#include "graph.h"
#include <cmath>

// Class representing a graph
Graph::Graph(ParityCheckMatrix& matrix) : matrix(matrix){
    std::vector<std::vector<int>> pcm_matrix = matrix.getBinaryMatrix();
    int i = 0;
    equalityNodesSize = pcm_matrix[i].size();
    checkNodesSize = pcm_matrix.size();
    adjListCheckNodes = std::vector<std::vector<std::pair<int, double>>>(checkNodesSize);
    adjListEqualityNodes = std::vector<std::vector<std::pair<int, double>>>(equalityNodesSize);
}

int Graph::getEqualityNodesSize(){
    return equalityNodesSize;
}

int Graph::sign(double x) {
    if (x > 0) {
        return 1;  // Positive
    } else if (x < 0) {
        return -1; // Negative
    } else {
        return 0;  // Zero
    }
}

void Graph::addEdge(int src, int dest)
{
    // Add a link from the check node to the equality node
    // At the beginning the vector inside the vector is empty, so we don't select any position of the vector
    // With push_back we add a new element to the end of the vector
    adjListCheckNodes[src].push_back({dest, 0});
    // Add a link from the equality node to the check node
    adjListEqualityNodes[dest].push_back({src, 0});
}

void Graph::generateGraph()
{
    std::vector<std::vector<int>> pcm_matrix = matrix.getBinaryMatrix();
    std::cout << "1" << std::endl;
    for (int i = 0; i < pcm_matrix.size(); ++i) {
        for (int j = 0; j < pcm_matrix[i].size(); ++j) {
            if (pcm_matrix[i][j] == 1) {
                addEdge(i, j); 
            }
        }
    }
}

double Graph::phi_tilde(double x) {
    // Calculate e^x
    double exp_x = std::exp(x);
    // Calculate the fraction (e^x + 1) / (e^x - 1)
    double fraction = (exp_x + 1) / (exp_x - 1);
    // Calculate the natural logarithm of the fraction
    return std::log(fraction);
}

std::vector<int> Graph::messagePassing(std::vector<double> receivedFromChannel, double variance){
    /*Messages that leave the equality nodes*/
    std::vector<double> vectorEqualityNodes(equalityNodesSize);
    /*Messages that leave the g nodes, coming from the channel*/
    std::vector<double> vectorgNodes(equalityNodesSize);
    double LLR_g;
    bool codeword = false;
    int counter = 0;
    std::vector<int> decodedBits(equalityNodesSize);

    for(int i = 0; i < receivedFromChannel.size(); ++i)
    {
        LLR_g = -(2.0 * receivedFromChannel[i]) / (variance * variance);
        vectorgNodes[i] = LLR_g;
    }

    while(!codeword){
        /*Messages that leave the equality nodes*/
        for(int i = 0; i < adjListEqualityNodes.size(); i++)
        {
            /*Check all the links connected to the current equality node*/
            for(int j = 0; j < adjListEqualityNodes[i].size(); ++j)
            {   
                /*Get the destination of the current link, so the check node*/
                int dest = adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the check nodes*/
                int sum = vectorgNodes[i];
                /*Check all the links connected to the current equality node*/
                for(int j_first = 0; j_first < adjListEqualityNodes[i].size(); j_first++){
                    /*If the current link is the same as the one we are checking, skip it*/
                    if(j_first == j){
                         continue;
                    }
                    int dest_first = adjListEqualityNodes[i][j_first].first;
                    /*Check all the links connected to the current check node*/
                    for(int j_second = 0; j_second < adjListCheckNodes[dest_first].size(); j_second++){
                        /*If destination of the current check nodes is the current equality node, add the message to the sum*/
                        if(adjListCheckNodes[dest_first][j_second].first == i){
                            sum += adjListCheckNodes[dest_first][j_second].second;
                            break;
                        }
                    }
                }
                /*Set the message to the current equality node*/
                adjListEqualityNodes[i][j].second = sum;
            }
        }

        /*Messages that leave the check nodes*/
        for(int i = 0; i < adjListCheckNodes.size(); i++)
        {
            /*Check all the links connected to the current check node*/
            for(int j = 0; j < adjListCheckNodes[i].size(); j++)
            {
                /*Get the destination, so the equality node*/
                int dest = adjListCheckNodes[i][j].first;
                /*Sum of the messages coming from the equality nodes except the destination*/
                double sum_of_LLR = 0;
                double product_sign = 1;
                for(int j_first = 0; j_first < adjListCheckNodes[i].size(); j_first++)
                {
                    if(j_first == j)
                    {
                        continue;
                    }
                    /*Get the destination, so the equality node*/
                    int dest_first = adjListCheckNodes[i][j_first].first;
                    for(int j_second = 0; j_second < adjListEqualityNodes[dest_first].size(); j_second++)
                    {
                        if(adjListEqualityNodes[dest_first][j_second].first == i)
                        {
                            sum_of_LLR += phi_tilde(std::fabs(adjListEqualityNodes[dest_first][j_second].second));
                            product_sign *= sign(adjListEqualityNodes[dest_first][j_second].second);
                            break;
                        }
                    }
                }
                /*Set the message to the current check node*/
                adjListCheckNodes[i][j].second = product_sign * phi_tilde(sum_of_LLR);
            }
        }

        /*Sum of the messages that enter in the equality nodes*/
        for(int i = 0; i < adjListEqualityNodes.size(); i++){
            double sum = 0;
            for(int j = 0; j < adjListEqualityNodes[i].size(); j++){
                /*Get the destination, so the check node*/
                int dest = adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the check nodes*/
                for(int j_first = 0; j_first < adjListCheckNodes[dest].size(); j_first++){
                    if(adjListCheckNodes[dest][j_first].first == i){
                        sum += adjListCheckNodes[dest][j_first].second;
                        break;
                    }
                }      
            }
            vectorEqualityNodes[i] = sum;
        }
    
        /*Marginalization*/
        for(int i = 0; i < vectorEqualityNodes.size(); i++){
            if(vectorEqualityNodes[i] + vectorgNodes[i] > 0){
                decodedBits[i] = 0;
            }else{
                decodedBits[i] = 1;
            }
        }
        if(matrix.isCodewordVector(decodedBits)){
            codeword = true;
        }
        if(counter > 100){
            break;
        }
    }
    return decodedBits;
}

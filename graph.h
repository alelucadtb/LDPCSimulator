// C++ Program to Implement Adjacency List
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include "paritycheckmatrix.h"

// Class representing a graph
class Graph {
    
public:
    
    /*Function to calculate the phi_tilde function*/
    double phi_tilde(double x);
    /*Function to calculate the sign of a number*/
    double sign(double x);

    // Constructor to initialize the graph
    // Parameters: vertices - number of vertices in the
    // graph
    //  directed - flag to indicate if the graph is directed
    //  (default is false)
    Graph(ParityCheckMatrix& pcm);
    
    /*Vector of vectors of pair vectors that represents the message sent from the equality nodes to the check nodes*/ 
    std::vector<std::vector<std::pair<int, double>>> adjListCheckNodes;
    /*Vector of vectors of pair vectors that represents the message sent from the check nodes to the equality nodes */ 
    std::vector<std::vector<std::pair<int, double>>> adjListEqualityNodes;
    /*Vector of vectors of pair vectors that represents the message sent from the w nodes to the equality nodes*/ 
    std::vector<std::vector<std::pair<int, double>>> adjListWNodes;
    
    ParityCheckMatrix& matrix;
    /*Number of check nodes*/
    int checkNodesSize;
    /*Number of equality nodes*/
    int equalityNodesSize;
    /*Get the number of equality nodes*/
    int getEqualityNodesSize();

    // Function to add an edge to the graph
    // Parameters: src - source vertex
    // dest - destination vertex
    void addEdge(int src, int dest);

    /*Generate the graph from the parity check matrix*/
    void generateGraph();

    /*Print the graph*/
    void printGraph(); 

    /*Execute the message passing algorithm on the graph*/
    std::vector<int> messagePassing(std::vector<double> receivedFromChannel, double variance);

    std::vector<int> upperMessagePassing(std::vector<double> LLRValues);
};

#endif
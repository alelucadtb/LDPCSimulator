#ifndef DECODER_H
#define DECODER_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <utility>
#include "paritycheckmatrix.h"
#include "graph.h"
#include "pam.h"

class Decoder{
    public:
        /*Constructor*/
        Decoder(std::vector<double> received_word, Graph& graph, double variance, std::vector<int> modulated_word, PAM& pam);
        // Constructor for the interleaving path
        Decoder(std::vector<double> received_word, Graph& graph, std::vector<double> varianceVector, std::vector<int> modulated_word, PAM& pam);

        /*Destructor*/
        ~Decoder();

        /*Decode the received word*/
        std::pair<std::vector<int>, int> BICMDecodingCycle(int fast_decoding_cycle); 
        /*Print the graph*/
        void printGraph();
        // For the interleaving path
        std::pair<std::vector<int>, int> interleavingBICMDecodingCycle(int fast_decoding_cycle); 

    private:
        /*The received word from the channel (r_l)*/
        std::vector<double> received_word;
        /*The graph related to the parity check matrix*/
        Graph& graph;
        /*The variance*/
        double variance;
        std::vector<double> varianceVector;
        /*The modulated word from the M-PAM modulator (d_l). It is what we send to the channel*/
        std::vector<int> modulated_word;
        /**
         * Function to add a link to the graph
         * @param src - source vertex
         * @param dest - destination vertex
        */
        void addLink(int src, int dest);
        /**
         * Generate the graph
        */
        void generateGraph();
        // Print the graph
        void printGraphBICM();
        /*Vector of vectors of pair vectors that represents the message sent from the w to the equality nodes*/   
        std::vector<std::vector<std::pair<int, double>>> adjListWNodes;
        /*Vector of vectors of pair vectors that represents the message sent from the equality nodes to the w nodes*/
        std::vector<std::pair<int, double>> adjListEqqNodes;
        /*Number of w nodes*/
        int wNodesSize;
        /*Number of equality nodes*/
        int equalityNodesSize;
        /*The PAM related to the modulated word*/
        PAM& pam;
        /*LLR for the w nodes*/
        std::vector<double> LLRwNodes;
        /*map that associates to a groups of bit the corresponding gray number*/
        int mapGrayNumber(std::vector<int> bit);
        /*LLR for the w nodes*/
        std::vector<double> calculateLLRwNodes(int output_link, int node_number);
        /*Calculate g function*/
        double calculateGFunction(int d, int i);
        // For the fact that we have different variance if we consider the Markov model
        double calculateInterleavingGFunction(int d, int i);
        /*Generate the permutation of log2(M)-1 bits*/
        std::vector<std::vector<int>> generatePermutation();
        /*Flip a bit*/
        int flipBit(int bit);
       
};

#endif
#ifndef DECODER_H
#define DECODER_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include "paritycheckmatrix.h"
#include "graph.h"
#include "pam.h"

class Decoder{
    public:
        /*Constructor*/
        Decoder(std::vector<double> received_word, Graph& graph, double variance, std::vector<int> modulated_word, PAM& pam);

        /*Destructor*/
        ~Decoder();

        /*Decode the received word*/
        std::vector<int> fastDecodingCycle();
    
    private:
        /*The received word from the channel (r_l)*/
        std::vector<double> received_word;
        /*The graph related to the parity check matrix*/
        Graph& graph;
        /*The variance*/
        double variance;
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
        /*Vector of vectors of pair vectors that represents the message sent from the w to the equality nodes*/   
        std::vector<std::vector<std::pair<int, double>>> adjListWNodes;
        /*Vector of vectors of pair vectors that represents the message sent from the equality nodes to the g nodes*/
        std::vector<std::vector<std::pair<int, double>>> adjListEqqNodes;
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
        std::vector<double> calculateLLRwNodes(int output_link);
        /*Calculate g function*/
        std::vector<double> calculateGFunction(int d, int i);
        /*Generate the permutation of log2(M)-1 bits*/
        std::vector<std::vector<int>> generatePermutation();
        /*Flip a bit*/
        int flipBit(int bit);
};

#endif
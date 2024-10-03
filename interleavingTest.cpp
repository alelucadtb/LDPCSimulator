#include <iostream>
#include "paritycheckmatrix.h"
#include "encoder.h"
#include "word.h"
#include "pam.h"
#include "channel.h"
#include "graph.h"
#include "decoder.h"
#include "error.h"
#include "interleaving.cpp"

std::vector<int> generateRandomVector(int len) {
    std::vector<int> random_vector(len);
    for (int i = 0; i < len; ++i) {
        random_vector[i] = rand() % 2;
    }
    return random_vector;
}

int main() {
    std::vector<double> iValue;
    int depths = 0;

    
    // Temporary value for the BER and SNR
    double BERtmp = 0.0;
    double SNRtmp = 0.0;
    int counter = 0;
    double R = 5.0/6.0;
    int count = 0;
    /*Create the parity check matrix*/
    ParityCheckMatrix pcm = ParityCheckMatrix(648, 27, R, "n648_Z27_R56.txt");
    pcm.writeBinaryMatrix("binary_matrix.txt");
    Encoder encoder = Encoder(pcm);
    std::vector<std::vector<int>> G = encoder.getGeneratingMatrix();
    std::vector<std::vector<int>> interleavingMatrix;
    
    for(int i = 0; i < depths; i++){
        /*Generate a random data word*/
        Word data_word = Word(generateRandomVector(540));
        /*Encode the data word*/
        Word encoded_word = encoder.encode(data_word);
    
        PAM modulated_pam = PAM(encoded_word, 8);
        /*Check if the encoded word is a codeword*/
        std::vector<int> modulated_word;

        if(pcm.isCodeword(encoded_word)){
            /*PAM Modulation*/
            modulated_word = modulated_pam.MPAMModulate(encoded_word);
            // This matrix has L number of rows that are the modulated word
            interleavingMatrix.push_back(modulated_word);
        }
    
    }
    /*for(int i = 0; i < modulated_word.size(); i++){
        std::cout << "modulated_world["<<i<<"]: " << modulated_word[i] << std::endl;
    }*/
    // Interleaving
    Interleaving inter = Interleaving();
    // Matrix for the interleaving
    std::vector<std::vector<double>> interleaved = inter.interleaving(interleavingMatrix);
    
    std::vector<std::vector<double>> received_block; 
    
    /*AWGN Channel*/
    Channel channel = Channel();
    for(int i = 0; i < interleaved.size(); i++){   
        std::vector<double> received_word = channel.markovChannel(interleaved[i]);
        received_block.push_back(received_word);
    }
    
    // Deinterleaving
    inter.deinterleaving(received_block);

    /*for(int i = 0; i < received_word.size(); i++){
        std::cout << "received_word["<< i <<"]: " << received_word[i] << std::endl;
    }*/
    
    /*Create the graph*/
    Graph graph = Graph(pcm);
    graph.generateGraph();
    //graph.printGraph();
    /*Execute the message passing algorithm*/
    //std::vector<int> decoded_word_2 = graph.messagePassing(received_word, 0.5);
    /*Print the decoded word*/
    Decoder decoder = Decoder(received_word, graph, numberOfCodeword[j].first, modulated_word, modulated_pam);
    std::vector<int> decoded_word_2 = decoder.BICMDecodingCycle(0);
    //std::cout << encoded_word << std::endl;
    
    //std::vector<int> decoded_word_2 = decoder.testingMethod(1);
    /*for(int i = 0; i < decoded_word_2.size(); i++){
        std::cout << decoded_word_2[i];
    }
    std::cout << std::endl;*/
    Error error = Error();
    BERtmp += error.calculateError(decoded_word_2, encoded_word);
    counter++;
    SNRtmp = (5.0/3.0) * 1.0/(numberOfCodeword[j].first * 2.0);
    //std::cout << SNRtmp << std::endl;
    std::cout << BERtmp/(double)numberOfCodeword[j].second << "\t" << SNRtmp << "\t" << numberOfCodeword[j].first << std::endl;
    std::cout << "counter: " << counter << std::endl;

    return 0;
}

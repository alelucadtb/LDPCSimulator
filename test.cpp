#include <iostream>
#include "paritycheckmatrix.h"
#include "encoder.h"
#include "word.h"
#include "pam.h"
#include "channel.h"
#include "graph.h"
#include "decoder.h"
#include "error.h"

std::vector<int> generateRandomVector(int len) {
    std::vector<int> random_vector(len);
    for (int i = 0; i < len; ++i) {
        random_vector[i] = rand() % 2;
    }
    return random_vector;
}

int main() {
    double R = 5.0/6.0;
    /*Create the parity check matrix*/
    ParityCheckMatrix pcm = ParityCheckMatrix(648, 27, R, "n648_Z27_R56.txt");
    pcm.writeBinaryMatrix("binary_matrix.txt");

    Encoder encoder = Encoder(pcm);
    std::vector<std::vector<int>> G = encoder.getGeneratingMatrix();
    
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
    }

    //for(int i = 0; i < modulated_word.size(); i++){
        //std::cout << "modulated_world["<<i<<"]: " << modulated_word[i] << std::endl;
    //}
    
    /*AWGN Channel*/
    AWGN awgn = AWGN(0.0, 0.5, modulated_word.size());
    Channel channel = Channel(awgn);
    std::vector<double> received_word = channel.AWGNChannel(modulated_word);
    /*for(int i = 0; i < received_word.size(); i++){
        std::cout << "received_word["<< i <<"]: " << received_word[i] << std::endl;
    }*/
    
    /*Create the graph*/
    Graph graph = Graph(pcm);
    graph.generateGraph();
    //graph.printGraph();
    /*Execute the message passing algorithm*/
    //std::vector<int> decoded_word_2 = graph.messagePassing(received_word, 0.25);
    /*Print the decoded word*/
    Decoder decoder = Decoder(received_word, graph, 0.5, modulated_word, modulated_pam);
    std::vector<int> decoded_word_2 = decoder.BICMDecodingCycle(0);
    std::cout << encoded_word << std::endl;
    
    //std::vector<int> decoded_word_2 = decoder.testingMethod(0.25);
    for(int i = 0; i < decoded_word_2.size(); i++){
        std::cout << decoded_word_2[i];
    }
    std::cout << std::endl;

    Error error = Error();
    std::cout << error.calculateError(decoded_word_2, encoded_word) << std::endl;
    return 0;
}   
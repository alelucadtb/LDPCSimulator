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
    std::vector<double> iValue;
    // Associate for different value of the variance a different number of repetion of the cycle
    std::vector<std::pair<double, int>> numberOfCodeword = {{0.3, 10000}, {0.4, 10000}, {0.5, 10000}, {0.6, 5000}, {0.7, 1000}, {0.8, 1000}, {1.0, 1000}, {1.1, 100}, {1.2, 10}, {1.3, 10}, {1.4, 10}, {1.5, 10}, {2, 10}, {2.5, 10}, {3, 10}, {3.5, 10}};
    //std::vector<std::pair<double, int>> numberOfCodeword = {{0.3, 1}, {0.4, 1}, {0.5, 1}, {0.51, 1}, {0.53, 1}, {0.54, 1}, {0.56, 1}, {0.57, 1}, {0.6, 1}, {0.7, 1}, {0.8, 1}};

    for(int j = 0; j < numberOfCodeword.size(); j++){
        // Temporary value for the BER and SNR
        double BERtmp = 0.0;
        double SNRtmp = 0.0;
        int counter = 0;
        for(int z = 0; z < numberOfCodeword[j].second; z++){
            double R = 5.0/6.0;
            int count = 0;
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

            /*for(int i = 0; i < modulated_word.size(); i++){
                std::cout << "modulated_world["<<i<<"]: " << modulated_word[i] << std::endl;
            }*/
            
            /*AWGN Channel*/
            AWGN awgn = AWGN(0.0, numberOfCodeword[j].first, modulated_word.size());
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
        }
        SNRtmp = (5.0/3.0) * 1.0/(numberOfCodeword[j].first * 2.0);
        //std::cout << SNRtmp << std::endl;
        std::cout << BERtmp/(double)numberOfCodeword[j].second << "\t" << SNRtmp << "\t" << numberOfCodeword[j].first << std::endl;
        std::cout << "counter: " << counter << std::endl;
    }

    return 0;
}

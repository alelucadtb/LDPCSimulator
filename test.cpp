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
    std::vector<double> BER;
    std::vector<double> SNR;
    int counter = 0;
    for(int j = 0; j < 10; j++){
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
        
        for(double i = 0.25; i < 3.5; i = i + 0.1){
            /*AWGN Channel*/
            AWGN awgn = AWGN(0.0, i, modulated_word.size());
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
            Decoder decoder = Decoder(received_word, graph, i, modulated_word, modulated_pam);
            std::vector<int> decoded_word_2 = decoder.BICMDecodingCycle(1);
            //std::cout << encoded_word << std::endl;
            
            //std::vector<int> decoded_word_2 = decoder.testingMethod(1);
            /*for(int i = 0; i < decoded_word_2.size(); i++){
                std::cout << decoded_word_2[i];
            }
            std::cout << std::endl;*/

            Error error = Error();
            if(counter == 0){
                BER.push_back(error.calculateError(decoded_word_2, encoded_word));
                SNR.push_back(1/(i*2));
            }else{
                BER[count] += error.calculateError(decoded_word_2, encoded_word);
            }
            std::cout << "count: " << count << std::endl;
            count++;
        }
        std::cout << "counter: " << counter << std::endl;
        counter++;
    }
    for(int i = 0; i < BER.size(); i++){
        std::cout <<  BER[i]/counter << "\t" << SNR[i] << std::endl;
    }
    std::cout << "counter: " << counter << std::endl;
    return 0;
}

#include <iostream>
#include "paritycheckmatrix.h"
#include "encoder.h"
#include "word.h"
#include "pam.h"
#include "channel.h"
#include "graph.h"
#include "decoder.h"
#include "error.h"
#include "interleaving.h"

std::vector<int> generateRandomVector(int len) {
    std::vector<int> random_vector(len);
    for (int i = 0; i < len; ++i) {
        random_vector[i] = rand() % 2;
    }
    return random_vector;
}

std::vector<std::vector<double>> divideVector(const std::vector<double>& original, size_t size) {
    std::vector<std::vector<double>> results;
    
    for (size_t i = 0; i < original.size(); i += size) {
        // Create a new subvector
        std::vector<double> subvector(original.begin() + i, original.begin() + std::min(original.size(), i + size));
        results.push_back(subvector);
    }
    
    return results;
}

int main() {
    std::vector<int> depthVector = {1, 5, 10, 20, 50, 100, 200};
    
    for(int i = 0; i < depthVector.size(); i++){
        int depths = depthVector[i];
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
        std::vector<std::vector<std::vector<int>>> tmpInterleavingMatrix(depths);
        // These vectors are useful for saving the possible words
        std::vector<std::vector<PAM>> modulated_pam_vector(depths);
        std::vector<std::vector<std::vector<int>>> modulated_word_vector(depths);
        std::vector<std::vector<Word>> encoded_word_vector(depths);
        int codewordsInABlock = 200/depths;
        std::vector<int> blockOfWords;
        int position = 0;

        for(int i = 0; i < depths; i++){
            for(int j = 0; j < codewordsInABlock; j++){
                Word data_word = Word(generateRandomVector(540));
                /*Encode the data word*/
                Word encoded_word = encoder.encode(data_word);
                encoded_word_vector[position].push_back(encoded_word);
                PAM modulated_pam = PAM(encoded_word, 8);
                modulated_pam_vector[position].push_back(modulated_pam);
                std::vector<int> modulated_word;

                /*Check if the encoded word is a codeword*/
                if(pcm.isCodeword(encoded_word)){
                    /*PAM Modulation*/
                    modulated_word = modulated_pam.MPAMModulate(encoded_word);
                    modulated_word_vector[position].push_back(modulated_word);
                    // This matrix has L number of rows that are the modulated word 
                    tmpInterleavingMatrix[position].push_back(modulated_word);
                }
            }
            position++;
            /*Generate a random data word*/
        }
        
        /*for(int l = 0; l < interleavingMatrix.size(); l++){
            for(int z = 0; z < interleavingMatrix[l].size(); z++){
                std::cout << "interleavingMatrix[" << l << "][" << z << "]: " << interleavingMatrix[l][z] << std::endl;
            }
        }*/
        Channel channel = Channel();
        std::vector<std::vector<int>> finalResults;
        int mod_count = 0;
        for(int l = 0; l < tmpInterleavingMatrix[0].size(); l++){
            std::vector<std::vector<int>> interleavingMatrix(depths);
            for(int n = 0; n < tmpInterleavingMatrix.size(); n++){
                interleavingMatrix[n] = tmpInterleavingMatrix[n][l];
            }
            // Interleaving
            Interleaving inter = Interleaving();
            // Matrix for the interleaving
            std::vector<std::vector<double>> interleaved = inter.interleaving(interleavingMatrix);
            std::vector<std::vector<double>> received_block;
            // The chosen variance of the Markov channel, depending on the good state or bad state
            /*AWGN Channel*/
            
            for(int k = 0; k < interleaved.size(); k++){ 
                auto result = channel.markovChannel(interleaved[k]); // Call once
                std::vector<double> received_word = result.first; // Get the received word
                std::vector<double> tmpVariance = result.second; // Get the variance
                received_block.push_back(received_word);
                /*for(int h = 0; h < received_word.size(); h++){
                    std::cout << "Channel word: " << received_word[h] << " ";
                }
                std::cout << std::endl;
                for(int l = 0; l < received_word.size(); l++){
                    std::cout << "Channel variance: " << tmpVariance[l] << " ";
                }
                std::cout << std::endl;*/
            }

            std::vector<std::vector<double>> deinterleaved = inter.deinterleaving(received_block);
            /*Create the graph*/
            Graph graph = Graph(pcm);
            graph.generateGraph();

            for(int j = 0; j < deinterleaved.size(); j++){
                Decoder decoder = Decoder(deinterleaved[j], graph, 0.85, modulated_word_vector[j][mod_count], modulated_pam_vector[j][mod_count]);
                std::pair<std::vector<int>, int> decoded_word = decoder.BICMDecodingCycle(1);
                finalResults.push_back(decoded_word.first);
                
            }
            Error error = Error();
            // For every possible results done by the interleaving depth
            for(int i = 0; i < finalResults.size(); i++){
                BERtmp += error.calculateError(finalResults[i], encoded_word_vector[i][mod_count]);
                counter++;
            }
            finalResults.clear();
            mod_count++;
        }

        /*for(int z = 0; z < finalResults.size(); z++){
            for(int k = 0; k < finalResults[z].size(); k++){
                std::cout << finalResults[z][k] << " ";
            }
            std::cout << std::endl;
        }*/

        std::cout << "Total BER: " << BERtmp/counter << std::endl;
        std::cout << "Depth: " << depths << std::endl;
    }
    
    return 0;
}

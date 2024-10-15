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
        std::vector<double> subvector(original.begin() + i,
                                       original.begin() + std::min(original.size(), i + size));
        results.push_back(subvector);
    }
    
    return results;
}

int main() {
    std::vector<int> depthVector = {1, 5, 10, 50, 100, 300};
    
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
        std::vector<std::vector<int>> interleavingMatrix;
        // These vectors are useful for saving the possible words
        std::vector<PAM> modulated_pam_vector;
        std::vector<std::vector<int>> modulated_word_vector;
        std::vector<Word> encoded_word_vector;
        int codewordsInABlock = 300/depths;
        std::vector<int> blockOfWords;

        for(int i = 0; i < depths; i++){
            for(int j = 0; j < codewordsInABlock; j++){
                Word data_word = Word(generateRandomVector(540));
                /*Encode the data word*/
                Word encoded_word = encoder.encode(data_word);
                encoded_word_vector.push_back(encoded_word);
                PAM modulated_pam = PAM(encoded_word, 8);
                modulated_pam_vector.push_back(modulated_pam);
                std::vector<int> modulated_word;
                

                /*Check if the encoded word is a codeword*/
                if(pcm.isCodeword(encoded_word)){
                    /*PAM Modulation*/
                    modulated_word = modulated_pam.MPAMModulate(encoded_word);
                    modulated_word_vector.push_back(modulated_word);
                    // This matrix has L number of rows that are the modulated word
                    blockOfWords.insert(blockOfWords.end(), modulated_word.begin(), modulated_word.end());   
                }
            }
            interleavingMatrix.push_back(blockOfWords);
            blockOfWords.clear();
            /*Generate a random data word*/
        }
        
        // Interleaving
        Interleaving inter = Interleaving();
        // Matrix for the interleaving
        std::vector<std::vector<double>> interleaved = inter.interleaving(interleavingMatrix);
        std::vector<std::vector<double>> received_block;
        // The chosen variance of the Markov channel, depending on the good state or bad state
        std::vector<std::vector<double>> lineVariance;

        /*for(int i = 0; i < interleaved.size(); i++){
            for(int j = 0; j < interleaved[i].size(); j++){
                std::cout << interleaved[i][j] << std::endl;
            }
        }*/
        
        /*AWGN Channel*/
        Channel channel = Channel();
        for(int i = 0; i < interleaved.size(); i++){   
            auto result = channel.markovChannel(interleaved[i]); // Call once
            std::vector<double> received_word = result.first; // Get the received word
            std::vector<double> tmpVariance = result.second; // Get the variance
            lineVariance.push_back(tmpVariance);
            received_block.push_back(received_word);
        }
    
        // Deinterleaving
        std::vector<std::vector<double>> deinterleaved = inter.deinterleaving(received_block);
        std::vector<std::vector<double>> deLineVariance = inter.deinterleaving(lineVariance);
        
        /*Create the graph*/
        Graph graph = Graph(pcm);
        graph.generateGraph();
        
        // Vector for saving the results
        std::vector<std::vector<int>> finalResults;
        int mod_count = 0;

        /*Print the decoded word*/
        for(int i = 0; i < deinterleaved.size(); i++){
            std::vector<std::vector<double>> listOfCodewords = divideVector(deinterleaved[i], 216);
            std::vector<std::vector<double>> listOfVariances = divideVector(deLineVariance[i], 216);
            //std::cout << listOfCodewords.size() << " " << listOfVariances.size() << std::endl;
            for(int j = 0; j < listOfCodewords.size(); j++){
                Decoder decoder = Decoder(listOfCodewords[j], graph, listOfVariances[j], modulated_word_vector[mod_count], modulated_pam_vector[mod_count]);
                std::pair<std::vector<int>, int> decoded_word = decoder.interleavingBICMDecodingCycle(1);
                finalResults.push_back(decoded_word.first);
                mod_count++;
            } 
        }
    
        Error error = Error();
        // For every possible results done by the interleaving depth
        for(int i = 0; i < finalResults.size(); i++){
            BERtmp += error.calculateError(finalResults[i], encoded_word_vector[i]);
            counter++;
        }
        std::cout << "Total BER: " << BERtmp/counter << std::endl;
        std::cout << "Depth: " << depths << std::endl;
    }
    
    return 0;
}

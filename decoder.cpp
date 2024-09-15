#include "decoder.h"

Decoder::Decoder(std::vector<double> received_word, Graph& graph, double variance, std::vector<int> modulated_word, PAM& pam) : received_word(received_word), graph(graph), variance(variance), modulated_word(modulated_word), pam(pam) {
    equalityNodesSize = graph.getEqualityNodesSize();
    wNodesSize =  received_word.size();
    adjListWNodes = std::vector<std::vector<std::pair<int, double>>>(wNodesSize);
    adjListEqqNodes = std::vector<std::vector<std::pair<int, double>>>(equalityNodesSize);
}

Decoder::~Decoder() {}

void Decoder::addLink(int src, int dest){
    adjListWNodes[src].push_back({dest, 0});
    adjListEqqNodes[dest].push_back({src, 0});
}

/* Generate the graph between the w-nodes and the equality nodes*/
void Decoder::generateGraph(){
    int bit_per_symbol = log2(pam.getM());
    for(int i = 0; i < wNodesSize; i = i + bit_per_symbol){
        for(int j = 0; j < bit_per_symbol; j++){
            /* Add a link between the w-node and the equality node */
            addLink(i, i+j);
        }
    }
}

int Decoder::flipBit(int bit){
    if(bit == 0){
        return 1;
    }
    return 0;
}

int Decoder::mapGrayNumber(std::vector<int> bit){
    std::string bitString;
    /*Convert the bit vector into a string*/
    for (int i = 0; i < bit.size(); ++i) {
        bitString += std::to_string(bit[i]);
    }
    Gray gray;
    int bit_per_symbol = log2(pam.getM());
    std::string gray_total = "";
    /*Generate the gray codes*/
    std::vector<std::string> gray_codes = gray.generateGrayarr(bit_per_symbol);
    /*Generate the gray codes in integer*/
    std::vector<int> gray_codes_int(gray_codes.size());
    for (int i = 0; i < gray_codes.size(); i++) {
        gray_codes_int[i] = 2*(i)-(pam.getM()-1);
    }
    /*Find the corresponding gray number*/
    for(int i = 0; i < gray_codes.size(); i++){
        if(gray_codes[i] == bitString){
            return gray_codes_int[i];
        }
    }
    return -1;
}

std::vector<double> Decoder::calculateLLRwNodes(int output_link){
    std::vector<double> LLRwNodes;
    for(int i = 0; i < adjListWNodes.size(); i++){
        for(int j = 0; j < adjListWNodes[i].size(); j++){
            int dest = adjListWNodes[i][j].first;
            for(int j_first = 0; j_first < adjListEqqNodes[dest].size(); j_first++){
                int dest_first = adjListEqqNodes[dest][j_first].first;
                if(dest_first == i || j_first != output_link){
                    LLRwNodes.push_back(adjListWNodes[i][j].second);
                }
            }
        }
    }
    return LLRwNodes;
}

std::vector<double> Decoder::calculateGFunction(int d, int l){
    std::vector<double> g;
    double temp_value_1 = std::fabs(received_word[l] - d);
    double temp_value_2 = variance*variance;
    g.push_back(exp(-temp_value_1*temp_value_1)/(2.0*(temp_value_2*temp_value_2)));
    return g;
}

std::vector<std::vector<int>> Decoder::generatePermutation(){
    int bit = log2(pam.getM()) - 1;
    std::vector<std::vector<int>> possible_permutation;
    std::vector<int> bit_permutation;
    std::string binary;
    for(int i = 0; i < bit; i++){
        binary = std::bitset<8>(i).to_string(); //to binary
        for (char& c : binary) {
            /*Convert the char into int*/
            bit_permutation.push_back(c - '0');
        }
        possible_permutation.push_back(bit_permutation);
        bit_permutation.clear();
    }
    return possible_permutation;
}

std::vector<int> Decoder::fastDecodingCycle(){
    std::vector<double> g;
    std::vector<double> LLRwNodes;
    std::vector<std::vector<int>> vector_map_0;
    std::vector<std::vector<int>> vector_map_1;
    std::vector<std::vector<int>> exponent_map_0;
    std::vector<std::vector<int>> exponent_map_1;
    double temp_1;
    std::vector<int> result_vector;
    //Generate the graph in the object of the class Decoder
    generateGraph();

    //Create a temporary vector_map element
    std::vector<int> vector_map_temp (log2(pam.getM()), 0);
    std::vector<int> exponent_map_temp (log2(pam.getM())-1, 0);
    //Create an array of possible permutation of M-1 bits
    std::vector<std::vector<int>> possible_permutation = generatePermutation();
    int bit_per_symbol = log2(pam.getM());
    //Generate the vector_map
    /* Select the position of the fixed bit in the vector_map */
    for(int i = 0; i < bit_per_symbol; i++){
        /* Select the value of the fixed bit in the vector_map */
        for(int i_first = 0; i_first < 2; i_first++){
            vector_map_temp[i] = i_first;
            for(int j = 0; j < possible_permutation.size(); j++){
                for(int j_first = 0; j_first < possible_permutation[j].size(); j_first++){
                    /* Flip the bit*/
                    exponent_map_temp[j_first] = flipBit(possible_permutation[j][j_first]);
                    if(i > j_first){
                        vector_map_temp[j_first] = possible_permutation[j][j_first];
                    }
                    /* If the position of the fixed bit in the vector_map is less than the position of the fixed bit in the possible_permutation */
                    else if(i < j_first){
                        vector_map_temp[j_first+1] = possible_permutation[j][j_first];
                    }
                }
                for(int j = 0; j < vector_map_temp.size(); j++){
                    std::cout << vector_map_temp[j] << " ";
                }
                std::cout << std::endl;
                if(i_first == 0){
                    vector_map_0.push_back(vector_map_temp);
                    exponent_map_0.push_back(exponent_map_temp);
                }
                else{
                    vector_map_1.push_back(vector_map_temp);
                    exponent_map_1.push_back(exponent_map_temp);
                }
                vector_map_temp.clear();
                exponent_map_temp.clear();
            }
            
        }
    }
    std::cout << "vector_map_0: " << vector_map_0.size() << std::endl;
    std::cout << "vector_map_1: " << vector_map_1.size() << std::endl;
   
    return result_vector;

}

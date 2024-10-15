#include "decoder.h"

/*double levels_100_5[] = {3.689088,2.996565,2.592140,2.305911,2.084631,1.904581,1.753105,1.622648,1.508333,1.406829,1.315757,1.233358,1.158296,1.089527,1.026224,0.967718,0.913459,0.862989,0.815924,0.771937,0.730746,0.692107,0.655809,0.621665,0.589509,0.559193,0.530588,0.503572,0.478041,0.453896,0.431048,0.409418,0.388929,0.369515,0.351111,0.333659,0.317105,0.301398,0.286492,0.272341,0.258906,0.246148,0.234031,0.222520,0.211584,0.201194,0.191320,0.181936,0.173017,0.164540,0.156482,0.148821,0.141538,0.134614,0.128030,0.121770,0.115818,0.110158,0.104775,0.099657,0.094789,0.090159,0.085757,0.081570,0.077587,0.073800,0.070198,0.066771,0.063513,0.060413,0.057465,0.054661,0.051994,0.049457,0.047044,0.044749,0.042566,0.040489,0.038514,0.036635,0.034848,0.033148,0.031531,0.029993,0.028530,0.027139,0.025815,0.024556,0.023358,0.022219,0.021135,0.020104,0.019124,0.018191,0.017304,0.016460,0.015657,0.014893,0.014167,0.013476};

double phi_tilde_quant(double x){

    int index = (int)(x / 0.05);
            if(index >= 100){
                return 0;
    }
    return levels_100_5[index];
}*/

Decoder::Decoder(std::vector<double> received_word, Graph& graph, double variance, std::vector<int> modulated_word, PAM& pam) : received_word(received_word), graph(graph), variance(variance), modulated_word(modulated_word), pam(pam) {
    equalityNodesSize = graph.getEqualityNodesSize();
    wNodesSize =  received_word.size();
    adjListWNodes = std::vector<std::vector<std::pair<int, double>>>(wNodesSize);
    
}

Decoder::Decoder(std::vector<double> received_word, Graph& graph, std::vector<double> varianceVector, std::vector<int> modulated_word, PAM& pam) : received_word(received_word), graph(graph), varianceVector(varianceVector), modulated_word(modulated_word), pam(pam) {
    equalityNodesSize = graph.getEqualityNodesSize();
    wNodesSize =  received_word.size();
    adjListWNodes = std::vector<std::vector<std::pair<int, double>>>(wNodesSize);
}

Decoder::~Decoder() {}

void Decoder::addLink(int src, int dest){
    adjListWNodes[src].push_back({dest, 0});
    adjListEqqNodes.push_back({src, 0});
}

/* Generate the graph between the w-nodes and the equality nodes*/
void Decoder::generateGraph(){
    int bit_per_symbol = log2(pam.getM());
    for(int i = 0; i < wNodesSize; i++){
        for(int j = 0; j < bit_per_symbol; j++){
            /* Add a link between the w-node and the equality node */
            addLink(i, i*bit_per_symbol + j);
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

std::vector<double> Decoder::calculateLLRwNodes(int output_link, int node_number){
    std::vector<double> LLRwNodes;
    // Checking all the link connected to the w nodes selected by the node_number
    for(int j = 0; j < adjListWNodes[node_number].size(); j++){
        int dest = adjListWNodes[node_number][j].first;
        // Checking all the link connected to the equality node selected by dest
        int dest_first = adjListEqqNodes[dest].first;
        // Checking if the equality node selected above has as a destination the w node selected by the node_number
        if(j != output_link){
            LLRwNodes.push_back(adjListEqqNodes[dest].second);
        }
    }
    /*for(int i = 0; i < LLRwNodes.size(); i++){
        std::cout << LLRwNodes[i] << " ";
    }
    std::cout << std::endl;*/
    return LLRwNodes;
}

void Decoder::printGraph(){
    for(int i = 0; i < adjListWNodes.size(); i++){
        for(int j = 0; j < adjListWNodes[i].size(); j++){
            std::cout << "adjListWNodes[" << i << "][" << j << "].first: " << adjListWNodes[i][j].first << std::endl;
            std::cout << "adjListWNodes[" << i << "][" << j << "].second: " << adjListWNodes[i][j].second << std::endl;
        }
    }
    for(int i = 0; i < adjListEqqNodes.size(); i++){
            std::cout << "adjListEqqNodes[" << i << "].first: " << adjListEqqNodes[i].first << std::endl;
            std::cout << "adjListEqqNodes[" << i << "].second: " << adjListEqqNodes[i].second << std::endl;
    }
}

double Decoder::calculateGFunction(int d, int l){
    double temp_value_1 = std::fabs(received_word[l] - d);
    return exp(-1.0 * pow(temp_value_1, 2.0)/(2.0 * variance));
}

double Decoder::calculateInterleavingGFunction(int d, int l){
    double temp_value_1 = std::fabs(received_word[l] - d);
    return exp(-1.0 * pow(temp_value_1, 2.0)/(2.0 * varianceVector[l]));
}

std::vector<std::vector<int>> Decoder::generatePermutation(){
    int bit = log2(pam.getM()) - 1;
    int num_options = pam.getM() / 2;
    std::vector<std::vector<int>> possible_permutation;
    std::vector<int> bit_permutation;
    std::string binary;
    for(int i = 0; i < num_options; i++){
        binary = std::bitset<8>(i).to_string(); //to binary
        for(int k = binary.length() - bit; k < binary.length(); k++){
            bit_permutation.push_back(binary[k] - '0');
        }
        possible_permutation.push_back(bit_permutation);
        bit_permutation.clear();
    }
    return possible_permutation;
}

std::pair<std::vector<int>, int> Decoder::BICMDecodingCycle(int fast_decoding_cycle){
    // Counters for the complexity calculation
    int addition_counter = 0; //1
    int multiplication_counter = 0; //2
    int division_counter = 0; //3
    int exp_counter = 0; //4
    int log_counter = 0; //4
    int comparison_counter = 0; //1
    int lookup_counter = 0; //1
    // End of the counters for the complexity calculation
    // Iterations of the while
    int counter = 0;
    // Identifies the output link of the conform nodes
    std::vector<int> output_link;
    // Identifies the possible input of the g function
    std::vector<std::vector<int>> vector_map_0;
    std::vector<std::vector<int>> vector_map_1;
    //Identifies the possible values for the exponential in the LLR formula for the conform nodes
    std::vector<std::vector<int>> exponent_map_0;
    std::vector<std::vector<int>> exponent_map_1;
    double temp_1;
    std::vector<int> result_vector;
    int bit_per_symbol = log2(pam.getM());
    //Generate the graph in the object of the class Decoder
    generateGraph();

    //Create a temporary vector_map element
    std::vector<int> vector_map_temp (bit_per_symbol, 0);
    std::vector<int> exponent_map_temp (bit_per_symbol-1, 0);
    
    // Create an array of possible permutation of M-1 bits
    std::vector<std::vector<int>> possible_permutation = generatePermutation();
    for(int i = 0; i < bit_per_symbol; i++){
        // Select the value of the fixed bit in the vector_map
        for(int i_first = 0; i_first < 2; i_first++){
            vector_map_temp[i] = i_first;
            for(int j = 0; j < possible_permutation.size(); j++){
                int k = 0;
                for(int j_first = 0; j_first < possible_permutation[j].size(); j_first++){
                    /* Flip the bit */
                    exponent_map_temp[j_first] = !possible_permutation[j][j_first];
                    //exponent_map_temp[j_first] = possible_permutation[j][j_first]? 1 : -1;
                    if(k == i){
                        k++;
                    }
                    vector_map_temp[k] =  possible_permutation[j][j_first];
                    k++;
                }
                if(i_first == 0){
                    vector_map_0.push_back(vector_map_temp);
                    exponent_map_0.push_back(exponent_map_temp);
                    output_link.push_back(i);
                }
                else{
                    vector_map_1.push_back(vector_map_temp);
                    exponent_map_1.push_back(exponent_map_temp);
                }
            }
        }
    }
    /*std::cout << "vectormap0" << std::endl;
    for(int i = 0; i < vector_map_0.size(); i++){
        for(int j = 0; j < vector_map_0[i].size(); j++){
            std::cout << vector_map_0[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << "vectormap1" << std::endl;
    for(int i = 0; i < vector_map_1.size(); i++){
        for(int j = 0; j < vector_map_1[i].size(); j++){
            std::cout << vector_map_1[i][j];
        }
        std::cout << std::endl;
    }*/
    //graph.printGraph();

    double LLR_g;
    bool codeword = false;
    std::vector<int> decodedBits(graph.equalityNodesSize);
    double sum_1 = 0.0;
    // Vector that contains the sum of the messages that enter in the equality nodes for the marginalization
    std::vector<double> vectorEqualityNodes(graph.equalityNodesSize);
   

    while(!codeword){
        if((fast_decoding_cycle == 1 && counter == 0) || fast_decoding_cycle == 0){
            double num_1;
            double den_1;
            double num_2;
            double den_2;
            double num_3;
            double den_3;
            double num;
            double den;
            int count = 0;
            // For checking which the link is the output of the conform nodes
            int out = output_link[0];
            std::vector<std::pair<double, double>> temp_vector;

            // Calculate the LLR for the conform nodes
            for(int i = 0; i < received_word.size(); i++){
                addition_counter++;
                // Select one of the possibile combinations of vector_map
                for(int j = 0; j < vector_map_0.size(); j++){
                    addition_counter++;
                    // Select one of the possibile combinations of exponent_map
                    std::vector<double> LLRwNodes = calculateLLRwNodes(output_link[j], i);
                    /*if(counter > 0){
                        for(int i = 0; i < LLRwNodes.size(); i++){
                            std::cout << "LLRwNodes["<<i<<"]" << LLRwNodes[i] << std::endl;
                        }
                        
                    }*/
                    //std::cout << "step 1" << std::endl;
                    num_1 = calculateGFunction(mapGrayNumber(vector_map_0[j]), i);
                    exp_counter++;
                    multiplication_counter = multiplication_counter + 2;
                    division_counter++;
                    den_1 = calculateGFunction(mapGrayNumber(vector_map_1[j]), i);
                    exp_counter++;
                    multiplication_counter = multiplication_counter + 2;
                    division_counter++;
                    num_2 = den_2 = 0.0;
                    //std::cout << "step 2" << std::endl;
                    //std::cout << "LLRwNodes size: " << LLRwNodes.size() << std::endl;
                    //std::cout << "exponent_map_0[0]:" << exponent_map_0[1].size() << std::endl;
                    for(int k = 0; k < exponent_map_0[j].size(); k++){
                        addition_counter++;
                        num_2 += ((double)exponent_map_0[j][k]) * LLRwNodes[k];
                        addition_counter++;
                        den_2 += ((double)exponent_map_1[j][k]) * LLRwNodes[k];
                        addition_counter++; 
                    }
                    
                    num_3 = den_3 = 0.0;
                    num_3 = exp(num_2);
                    exp_counter++;
                    den_3 = exp(den_2);
                    exp_counter++;
                    
                    if(out != output_link[j]){
                        temp_vector.push_back(std::make_pair(num, den));
                        out = output_link[j];
                        num = num_1 * num_3;
                        multiplication_counter++;
                        den = den_1 * den_3;
                        multiplication_counter++;
                    } else{
                        num += num_1 * num_3;
                        addition_counter++;
                        multiplication_counter++;
                        den += den_1 * den_3;
                        addition_counter++;
                        multiplication_counter++;
                    }

                }
                temp_vector.push_back(std::make_pair(num, den));
                out = output_link[0];
                num = den = 0.0;
                //std::cout << "temp_vector.size(): " << temp_vector.size() << std::endl;
                /*Calculate the LLR for the w-nodes*/
                for(int j = 0; j < temp_vector.size(); j++){
                    addition_counter++;
                    double LLR_temp = temp_vector[j].first/temp_vector[j].second;
                    division_counter++;
                    //std::cout << "LLR_temp: " << LLR_temp << std::endl;
                    adjListWNodes[i][j].second = log(LLR_temp);
                    log_counter++;
                    //std::cout << "adjListWNodes["<<i<<"]["<<j<<"]: " << adjListWNodes[i][j].second << std::endl;
                }
                //std::cout << "adjListWNodes[i].size(): " << adjListWNodes[i].size() << std::endl;
                temp_vector.clear();
            }
        }

        /*Messages that leave the equality nodes*/
        for(int i = 0; i < graph.adjListEqualityNodes.size(); i++)
        {
            addition_counter++;
            /*Check all the links connected to the current equality node*/
            for(int j = 0; j < graph.adjListEqualityNodes[i].size(); ++j)
            {   
                addition_counter++;
                double sum = 0.0;
                /*Get the destination of the current link, so the check node*/
                int dest = graph.adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the w-nodes to the equality node*/
                for(int l = 0; l < adjListWNodes.size(); l++){
                    addition_counter++;
                    for(int z = 0; z < adjListWNodes[l].size(); z++){
                        addition_counter++;
                        comparison_counter++;
                        if(adjListWNodes[l][z].first == i){
                            sum_1 = adjListWNodes[l][z].second;
                            //std::cout << "sum_1: " << sum_1 << std::endl;
                            goto jump;
                        }
                    }
                }
                jump:
                /*Check all the links connected to the current equality node*/
                for(int j_first = 0; j_first < graph.adjListEqualityNodes[i].size(); j_first++){
                    addition_counter++;
                    /*If the current link is the same as the one we are checking, skip it*/
                    comparison_counter++;
                    if(j_first == j){
                        continue;
                    }
                    int dest_first = graph.adjListEqualityNodes[i][j_first].first;
                    /*Check all the links connected to the current check node*/
                    for(int j_second = 0; j_second < graph.adjListCheckNodes[dest_first].size(); j_second++){
                        addition_counter++;
                        /*If destination of the current check nodes is the current equality node, add the message to the sum*/
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest_first][j_second].first == i){
                            addition_counter++;
                            sum += graph.adjListCheckNodes[dest_first][j_second].second;
                            break;
                        }
                    }
                }
                //std::cout << "sum: " << sum << std::endl;
                /*Set the message to the current equality node*/
                graph.adjListEqualityNodes[i][j].second = sum + sum_1;
            }
        }
        
        /*Messages that leave the check nodes*/
        for(int i = 0; i < graph.adjListCheckNodes.size(); i++)
        {
            addition_counter++;
            /*Check all the links connected to the current check node*/
            for(int j = 0; j < graph.adjListCheckNodes[i].size(); j++)
            {
                addition_counter++;
                /*Sum of the messages coming from the equality nodes except the destination*/
                double sum_of_LLR = 0.0;
                double product_sign = 1.0;
                for(int j_first = 0; j_first < graph.adjListCheckNodes[i].size(); j_first++){
                    addition_counter++;
                    comparison_counter++;
                    if(j_first == j)
                    {
                        continue;
                    }
                    /*Get the destination, so the equality node*/
                    int dest_first = graph.adjListCheckNodes[i][j_first].first;
                    for(int j_second = 0; j_second < graph.adjListEqualityNodes[dest_first].size(); j_second++)
                    {
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListEqualityNodes[dest_first][j_second].first == i)
                        {   
                            double modul_LLR = std::fabs(graph.adjListEqualityNodes[dest_first][j_second].second);
                            sum_of_LLR += graph.phi_tilde(modul_LLR);
                            addition_counter++;
                            product_sign *= graph.sign(graph.adjListEqualityNodes[dest_first][j_second].second);
                            multiplication_counter++;
                            break;   
                        }
                    }
                }
                /*Set the message to the current check node*/
                graph.adjListCheckNodes[i][j].second = product_sign * graph.phi_tilde(sum_of_LLR);
                multiplication_counter++;
            }   
        }

        //exit(0);
        /*Sum of the messages that enter in the equality nodes*/
        comparison_counter++;
        if(fast_decoding_cycle == 0){
            for(int i = 0; i < graph.adjListEqualityNodes.size(); i++){
                addition_counter++;
                double sum_5 = adjListWNodes[i/bit_per_symbol][i%bit_per_symbol].second;
                division_counter = division_counter + 2;
                double sum_4 = 0.0;
                for(int j = 0; j < graph.adjListEqualityNodes[i].size(); j++){
                    addition_counter++;
                    /*Get the destination, so the check node*/
                    int dest = graph.adjListEqualityNodes[i][j].first;
                    /*Sum of the messages coming from the check nodes*/
                    for(int j_first = 0; j_first < graph.adjListCheckNodes[dest].size(); j_first++){
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest][j_first].first == i){
                            sum_4 += graph.adjListCheckNodes[dest][j_first].second;
                            addition_counter++;
                            break;
                        }
                    }        
                }
                adjListEqqNodes[i].second = sum_4;
                vectorEqualityNodes[i] = sum_4 + sum_5; 
                addition_counter++;
            }
        }
        else{
            for(int i = 0; i < graph.adjListEqualityNodes.size(); i++){
                addition_counter++;
                double sum_3 = adjListWNodes[i/bit_per_symbol][i%bit_per_symbol].second;
                division_counter = division_counter + 2;
                for(int j = 0; j < graph.adjListEqualityNodes[i].size(); j++){
                    addition_counter++;
                    /*Get the destination, so the check node*/
                    int dest = graph.adjListEqualityNodes[i][j].first;
                    /*Sum of the messages coming from the check nodes*/
                    for(int j_first = 0; j_first < graph.adjListCheckNodes[dest].size(); j_first++){
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest][j_first].first == i){
                            addition_counter++;
                            sum_3 += graph.adjListCheckNodes[dest][j_first].second;
                        }
                    }      
                }
                vectorEqualityNodes[i] = sum_3;
            }
        }

        /*for(int z = 0; z < vectorEqualityNodes.size(); z++){
            std::cout << "vectorEqualityNodes[" << z << "]: " << vectorEqualityNodes[z] << std::endl;
        }*/
        /*Marginalization*/
        for(int i = 0; i < vectorEqualityNodes.size(); i++){
            addition_counter++;
            comparison_counter++;
            if(vectorEqualityNodes[i] > 0){
                decodedBits[i] = 0;
            }else{
                decodedBits[i] = 1;
            }
        }
        if(graph.matrix.isCodewordVector(decodedBits)){
            codeword = true;
        }
        addition_counter++;
        counter++;
        //std::cout << "counter: " << counter << std::endl;
        comparison_counter++;
        if(counter == 1){
            break;
        }
    }
    graph.printGraph();
    printGraph();
    std::cout << adjListEqqNodes.size() << std::endl;

    // Part for the complexity of the decoding
    // The cost for all the operations
    int cost_of_operations = addition_counter + multiplication_counter*2 + division_counter*3 + exp_counter*4 + log_counter*4 + comparison_counter;
    int number_of_operations = addition_counter + multiplication_counter + division_counter + exp_counter + log_counter + comparison_counter;
    //std::cout << "medium cost for an operation: " << cost_of_operations << std::endl;
    
    return std::make_pair(decodedBits, cost_of_operations);
}

std::pair<std::vector<int>, int> Decoder::interleavingBICMDecodingCycle(int fast_decoding_cycle){
    // Counters for the complexity calculation
    int addition_counter = 0; //1
    int multiplication_counter = 0; //2
    int division_counter = 0; //3
    int exp_counter = 0; //4
    int log_counter = 0; //4
    int comparison_counter = 0; //1
    int lookup_counter = 0; //1
    // End of the counters for the complexity calculation
    // Iterations of the while
    int counter = 0;
    std::vector<double> g;
    std::vector<int> output_link;
    std::vector<std::vector<int>> vector_map_0;
    std::vector<std::vector<int>> vector_map_1;
    std::vector<std::vector<int>> exponent_map_0;
    std::vector<std::vector<int>> exponent_map_1;
    double temp_1;
    std::vector<int> result_vector;
    int bit_per_symbol = log2(pam.getM());
    //Generate the graph in the object of the class Decoder
    generateGraph();

    //Create a temporary vector_map element
    std::vector<int> vector_map_temp (bit_per_symbol, 0);
    std::vector<int> exponent_map_temp (bit_per_symbol-1, 0);
    
    // Create an array of possible permutation of M-1 bits
    std::vector<std::vector<int>> possible_permutation = generatePermutation();
    for(int i = 0; i < bit_per_symbol; i++){
        // Select the value of the fixed bit in the vector_map
        for(int i_first = 0; i_first < 2; i_first++){
            vector_map_temp[i] = i_first;
            for(int j = 0; j < possible_permutation.size(); j++){
                int k = 0;
                for(int j_first = 0; j_first < possible_permutation[j].size(); j_first++){
                    /* Flip the bit */
                    exponent_map_temp[j_first] = !possible_permutation[j][j_first];
                    //exponent_map_temp[j_first] = possible_permutation[j][j_first]? 1 : -1;
                    if(k == i){
                        k++;
                    }
                    vector_map_temp[k] =  possible_permutation[j][j_first];
                    k++;
                }
                if(i_first == 0){
                    vector_map_0.push_back(vector_map_temp);
                    exponent_map_0.push_back(exponent_map_temp);
                    output_link.push_back(i);
                }
                else{
                    vector_map_1.push_back(vector_map_temp);
                    exponent_map_1.push_back(exponent_map_temp);
                }
            }
        }
    }

    //graph.printGraph();
    double LLR_g;
    bool codeword = false;
    std::vector<int> decodedBits(graph.equalityNodesSize);
    double sum_1 = 0;
    /* Vector that contains the sum of the messages that enter in the equality nodes*/
    std::vector<double> vectorEqualityNodes(graph.equalityNodesSize);
   

    while(!codeword){
        if((fast_decoding_cycle == 1 && counter == 0) || fast_decoding_cycle == 0){
            //std::cout << "if" << std::endl;
            double num_1;
            double den_1;
            double num_2;
            double den_2;
            double num;
            double den;
            int count = 0;
            // For checking which the link is the output
            int out = output_link[0];
            std::vector<std::pair<double, double>> temp_vector;

            for(int i = 0; i < received_word.size(); i++){
                addition_counter++;
                //Select one of the possibile combinations of vector_map
                for(int j = 0; j < vector_map_0.size(); j++){
                    addition_counter++;
                    //Select one of the possibile combinations of exponent_map
                    std::vector<double> LLRwNodes = calculateLLRwNodes(output_link[j], i);
                    /*if(counter > 0){
                        for(int i = 0; i < LLRwNodes.size(); i++){
                            std::cout << "LLRwNodes["<<i<<"]" << LLRwNodes[i] << std::endl;
                        }
                        
                    }*/
                    //std::cout << "step 1" << std::endl;
                    num_1 = calculateInterleavingGFunction(mapGrayNumber(vector_map_0[j]), i);
                    exp_counter++;
                    multiplication_counter = multiplication_counter + 2;
                    division_counter++;
                    den_1 = calculateInterleavingGFunction(mapGrayNumber(vector_map_1[j]), i);
                    exp_counter++;
                    multiplication_counter = multiplication_counter + 2;
                    division_counter++;
                    num_2 = den_2 = 0.0;
                    //std::cout << "step 2" << std::endl;
                    //std::cout << "LLRwNodes size: " << LLRwNodes.size() << std::endl;
                    //std::cout << "exponent_map_0[0]:" << exponent_map_0[1].size() << std::endl;
                    for(int k = 0; k < exponent_map_0[j].size(); k++){
                        addition_counter++;
                        num_2 += ((double)exponent_map_0[j][k]) * LLRwNodes[k];
                        addition_counter++;
                        den_2 += ((double)exponent_map_1[j][k]) * LLRwNodes[k];
                        addition_counter++; 
                    }
                    
                    num_2 = exp(num_2);
                    exp_counter++;
                    den_2 = exp(den_2);
                    exp_counter++;
                    
                    if(out != output_link[j]){
                        temp_vector.push_back(std::make_pair(num, den));
                        out = output_link[j];
                        num = num_1 * num_2;
                        multiplication_counter++;
                        den = den_1 * den_2;
                        multiplication_counter++;
                    } else{
                        num += num_1 * num_2;
                        addition_counter++;
                        multiplication_counter++;
                        den += den_1 * den_2;
                        addition_counter++;
                        multiplication_counter++;
                    }

                }
                temp_vector.push_back(std::make_pair(num, den));
                out = output_link[0];
                num = den = 0.0;
                //std::cout << "temp_vector.size(): " << temp_vector.size() << std::endl;
                /*Calculate the LLR for the w-nodes*/
                for(int j = 0; j < temp_vector.size(); j++){
                    addition_counter++;
                    double LLR_temp = temp_vector[j].first/temp_vector[j].second;
                    division_counter++;
                    //std::cout << "LLR_temp: " << LLR_temp << std::endl;
                    adjListWNodes[i][j].second = log(LLR_temp);
                    log_counter++;
                    //std::cout << "adjListWNodes["<<i<<"]["<<j<<"]: " << adjListWNodes[i][j].second << std::endl;
                }
                //std::cout << "adjListWNodes[i].size(): " << adjListWNodes[i].size() << std::endl;
                temp_vector.clear();
            }
        }

        /*Messages that leave the equality nodes*/
        for(int i = 0; i < graph.adjListEqualityNodes.size(); i++)
        {
            addition_counter++;
            /*Check all the links connected to the current equality node*/
            for(int j = 0; j < graph.adjListEqualityNodes[i].size(); ++j)
            {   
                addition_counter++;
                double sum = 0;
                /*Get the destination of the current link, so the check node*/
                int dest = graph.adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the w-nodes to the equality node*/
                for(int l = 0; l < adjListWNodes.size(); l++){
                    addition_counter++;
                    for(int z = 0; z < adjListWNodes[l].size(); z++){
                        addition_counter++;
                        comparison_counter++;
                        if(adjListWNodes[l][z].first == i){
                            sum_1 = adjListWNodes[l][z].second;
                            //std::cout << "sum_1: " << sum_1 << std::endl;
                            goto jump;
                        }
                    }
                }
                jump:
                /*Check all the links connected to the current equality node*/
                for(int j_first = 0; j_first < graph.adjListEqualityNodes[i].size(); j_first++){
                    addition_counter++;
                    /*If the current link is the same as the one we are checking, skip it*/
                    comparison_counter++;
                    if(j_first == j){
                        continue;
                    }
                    int dest_first = graph.adjListEqualityNodes[i][j_first].first;
                    /*Check all the links connected to the current check node*/
                    for(int j_second = 0; j_second < graph.adjListCheckNodes[dest_first].size(); j_second++){
                        addition_counter++;
                        /*If destination of the current check nodes is the current equality node, add the message to the sum*/
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest_first][j_second].first == i){
                            addition_counter++;
                            sum += graph.adjListCheckNodes[dest_first][j_second].second;
                            break;
                        }
                    }
                }
                //std::cout << "sum: " << sum << std::endl;
                /*Set the message to the current equality node*/
                graph.adjListEqualityNodes[i][j].second = sum + sum_1;
            }
        }
        
        /*Messages that leave the check nodes*/
        for(int i = 0; i < graph.adjListCheckNodes.size(); i++)
        {
            addition_counter++;
            /*Check all the links connected to the current check node*/
            for(int j = 0; j < graph.adjListCheckNodes[i].size(); j++)
            {
                addition_counter++;
                /*Sum of the messages coming from the equality nodes except the destination*/
                double sum_of_LLR = 0.0;
                double product_sign = 1.0;
                for(int j_first = 0; j_first < graph.adjListCheckNodes[i].size(); j_first++)
                {
                    addition_counter++;
                    comparison_counter++;
                    if(j_first == j)
                    {
                        continue;
                    }
                    /*Get the destination, so the equality node*/
                    int dest_first = graph.adjListCheckNodes[i][j_first].first;
                    for(int j_second = 0; j_second < graph.adjListEqualityNodes[dest_first].size(); j_second++)
                    {
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListEqualityNodes[dest_first][j_second].first == i)
                        {   
                            double modul_LLR = std::fabs(graph.adjListEqualityNodes[dest_first][j_second].second);
                            sum_of_LLR += graph.phi_tilde(modul_LLR);
                            addition_counter++;
                            product_sign *= graph.sign(graph.adjListEqualityNodes[dest_first][j_second].second);
                            multiplication_counter++;
                            break;   
                        }
                    }
                }
                /*Set the message to the current check node*/
                graph.adjListCheckNodes[i][j].second = product_sign * graph.phi_tilde(sum_of_LLR);
                multiplication_counter++;
            }   
        }

        //exit(0);
        /*Sum of the messages that enter in the equality nodes*/
        comparison_counter++;
        if(fast_decoding_cycle == 0){
            for(int i = 0; i < graph.adjListEqualityNodes.size(); i++){
                addition_counter++;
                double sum_5 = adjListWNodes[i/bit_per_symbol][i%bit_per_symbol].second;
                division_counter = division_counter + 2;
                double sum_4 = 0;
                for(int j = 0; j < graph.adjListEqualityNodes[i].size(); j++){
                    addition_counter++;
                    /*Get the destination, so the check node*/
                    int dest = graph.adjListEqualityNodes[i][j].first;
                    /*Sum of the messages coming from the check nodes*/
                    for(int j_first = 0; j_first < graph.adjListCheckNodes[dest].size(); j_first++){
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest][j_first].first == i){
                            sum_4 += graph.adjListCheckNodes[dest][j_first].second;
                            addition_counter++;
                        }
                    }        
                }
                adjListEqqNodes[i].second = sum_4;
                vectorEqualityNodes[i] = sum_4 + sum_5; 
                addition_counter++;
            }
        }
        else{
            for(int i = 0; i < graph.adjListEqualityNodes.size(); i++){
                addition_counter++;
                double sum_3 = adjListWNodes[i/bit_per_symbol][i%bit_per_symbol].second;
                division_counter = division_counter + 2;
                for(int j = 0; j < graph.adjListEqualityNodes[i].size(); j++){
                    addition_counter++;
                    /*Get the destination, so the check node*/
                    int dest = graph.adjListEqualityNodes[i][j].first;
                    /*Sum of the messages coming from the check nodes*/
                    for(int j_first = 0; j_first < graph.adjListCheckNodes[dest].size(); j_first++){
                        addition_counter++;
                        comparison_counter++;
                        if(graph.adjListCheckNodes[dest][j_first].first == i){
                            addition_counter++;
                            sum_3 += graph.adjListCheckNodes[dest][j_first].second;
                        }
                    }      
                }
                vectorEqualityNodes[i] = sum_3;
            }
        }
        /*Marginalization*/
        for(int i = 0; i < vectorEqualityNodes.size(); i++){
            addition_counter++;
            comparison_counter++;
            if(vectorEqualityNodes[i] > 0){
                decodedBits[i] = 0;
            }else{
                decodedBits[i] = 1;
            }
        }
        if(graph.matrix.isCodewordVector(decodedBits)){
            codeword = true;
        }
        addition_counter++;
        counter++;
        //std::cout << "counter: " << counter << std::endl;
        comparison_counter++;
        if(counter == 20){
            break;
        }
    }

    // Part for the complexity of the decoding
    // The cost for all the operations
    int cost_of_operations = addition_counter + multiplication_counter*2 + division_counter*3 + exp_counter*4 + log_counter*4 + comparison_counter;
    int number_of_operations = addition_counter + multiplication_counter + division_counter + exp_counter + log_counter + comparison_counter;
    //std::cout << "medium cost for an operation: " << cost_of_operations << std::endl;
    
    return std::make_pair(decodedBits, cost_of_operations);
}

std::vector<int> Decoder::fastDecodingCycle(){
    std::vector<double> g;
    std::vector<double> LLRwNodes;
    std::vector<int> output_link;
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
                int k = 0;
                for(int j_first = 0; j_first < possible_permutation[j].size(); j_first++){
                    /* Flip the bit */
                    exponent_map_temp[j_first] = possible_permutation[j][j_first];
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
   
    double num_1;
    double den_1;
    double num_2;
    double den_2;
    double num;
    double den;
    int out = output_link[0];
    std::vector<std::pair<double, double>> temp_vector;

    //Select the received symbols
    for(int i = 0; i < received_word.size(); i++){
        //Select one of the possibile combinations of vector_map
        for(int j = 0; j < vector_map_0.size(); j++){
            //Select one of the possibile combinations of exponent_map
            std::vector<double> LLRwNodes = calculateLLRwNodes(output_link[j], i);
            num_1 = calculateGFunction(mapGrayNumber(vector_map_0[j]), i);
            den_1 = calculateGFunction(mapGrayNumber(vector_map_1[j]), i);
            num_2 = den_2 = 0.0;
            for(int k = 0; k < exponent_map_0[j].size(); k++){
                for(int l = 0; l < LLRwNodes.size(); l++){
                    num_2 += exponent_map_0[j][k] * LLRwNodes[l];
                    den_2 += exponent_map_1[j][k] * LLRwNodes[l];
                }
            }
            num_2 = exp(num_2);
            den_2 = exp(den_2);
            num += num_1 * num_2;
            den += den_1 * den_2;
            if(out != output_link[j]){
                out = output_link[j];
                temp_vector.push_back(std::make_pair(num, den));
                num = 0.0;
                den = 0.0;
            }
        }
        /*Calculate the LLR for the w-nodes*/
        for(int j = 0; j < temp_vector.size(); j++){
            double LLR_temp = temp_vector[j].first/temp_vector[j].second;
            adjListWNodes[i][j].second = log(LLR_temp);
            std::cout << "adjListWNodes[i][j].second: " << adjListWNodes[i][j].second << std::endl;
        }
        temp_vector.clear();
    }
    
    graph.generateGraph();
    double LLR_g;
    bool codeword = false;
    int counter = 0;
    std::vector<int> decodedBits(graph.equalityNodesSize);
    double sum_1 = 0;
    double sum = 0;
    /* Vector that contains the sum of the messages that enter in the equality nodes*/
    std::vector<double> vectorEqualityNodes(graph.equalityNodesSize);
    /*Messages that leave the equality nodes*/

    while(!codeword){
        /*Messages that leave the equality nodes*/
        for(int i = 0; i < graph.adjListEqualityNodes.size(); i++)
        {
            /*Check all the links connected to the current equality node*/
            for(int j = 0; j < graph.adjListEqualityNodes[i].size(); ++j)
            {   
                sum = 0;
                /*Get the destination of the current link, so the check node*/
                int dest = graph.adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the w-nodes to the equality node*/
                for(int l = 0; l < adjListWNodes.size(); l++){
                    for(int z = 0; z < adjListWNodes[l].size(); z++){
                        if(adjListWNodes[l][z].first == i){
                            sum_1 = adjListWNodes[l][z].second;
                            std::cout << "sum_1: " << sum_1 << std::endl;
                            break;
                        }
                    }
                }
                /*Check all the links connected to the current equality node*/
                for(int j_first = 0; j_first < graph.adjListEqualityNodes[i].size(); j_first++){
                    /*If the current link is the same as the one we are checking, skip it*/
                    if(j_first == j){
                        continue;
                    }
                    int dest_first = graph.adjListEqualityNodes[i][j_first].first;
                    /*Check all the links connected to the current check node*/
                    for(int j_second = 0; j_second < graph.adjListCheckNodes[dest_first].size(); j_second++){
                        /*If destination of the current check nodes is the current equality node, add the message to the sum*/
                        if(graph.adjListCheckNodes[dest_first][j_second].first == i){
                            sum += graph.adjListCheckNodes[dest_first][j_second].second;
                            std::cout << "sum: " << sum << std::endl;
                            break;
                        }
                    }
                }
                /*Set the message to the current equality node*/
                graph.adjListEqualityNodes[i][j].second = sum + sum_1;
                //std::cout << "value: " << graph.adjListEqualityNodes[i][j].second << std::endl;
            }
        }

        /*Messages that leave the check nodes*/
        for(int i = 0; i < graph.adjListCheckNodes.size(); i++)
        {
            /*Check all the links connected to the current check node*/
            for(int j = 0; j < graph.adjListCheckNodes[i].size(); j++)
            {
                /*Get the destination, so the equality node*/
                int dest = graph.adjListCheckNodes[i][j].first;
                /*Sum of the messages coming from the equality nodes except the destination*/
                double sum_of_LLR = 0;
                double product_sign = 1;
                for(int j_first = 0; j_first < graph.adjListCheckNodes[i].size(); j_first++)
                {
                    if(j_first == j)
                    {
                        continue;
                    }
                    /*Get the destination, so the equality node*/
                    int dest_first = graph.adjListCheckNodes[i][j_first].first;
                    for(int j_second = 0; j_second < graph.adjListEqualityNodes[dest_first].size(); j_second++)
                    {
                        if(graph.adjListEqualityNodes[dest_first][j_second].first == i)
                        {   
                            double modul_LLR = std::fabs(graph.adjListEqualityNodes[dest_first][j_second].second);
                            sum_of_LLR += graph.phi_tilde(modul_LLR);
                            product_sign *= graph.sign(graph.adjListEqualityNodes[dest_first][j_second].second);
                            //std::cout << "value: " << graph.adjListEqualityNodes[dest_first][j_second].second << std::endl;
                            break;   
                        }
                    }
                }
                /*Set the message to the current check node*/
                graph.adjListCheckNodes[i][j].second = product_sign * graph.phi_tilde(sum_of_LLR);
            }   
        }
        //exit(0);
        /*Sum of the messages that enter in the equality nodes*/
        for(int i = 0; i < graph.adjListEqualityNodes.size(); i++){
            double sum_3 = adjListWNodes[i/bit_per_symbol][i%bit_per_symbol].second;
            for(int j = 0; j < graph.adjListEqualityNodes[i].size(); j++){
                /*Get the destination, so the check node*/
                int dest = graph.adjListEqualityNodes[i][j].first;
                /*Sum of the messages coming from the check nodes*/
                for(int j_first = 0; j_first < graph.adjListCheckNodes[dest].size(); j_first++){
                    if(graph.adjListCheckNodes[dest][j_first].first == i){
                        sum_3 += graph.adjListCheckNodes[dest][j_first].second;
                    }
                }      
            }
            vectorEqualityNodes[i] = sum_3;
        }
        /*Marginalization*/
        for(int i = 0; i < vectorEqualityNodes.size(); i++){
            if(vectorEqualityNodes[i] > 0){
                decodedBits[i] = 0;
            }else{
                decodedBits[i] = 1;
            }
        }
        if(graph.matrix.isCodewordVector(decodedBits)){
            codeword = true;
        }
        counter++;
        if(counter > 100){
            break;
        }
    }

    std::cout << "counter: " << counter << std::endl;
    return decodedBits;

}
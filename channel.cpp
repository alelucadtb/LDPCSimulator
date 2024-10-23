#include "channel.h"

Channel::Channel(AWGN& awgn) : awgn(awgn) {
}

Channel::Channel() {
    this->awgn = AWGN();
    this->actualState = 0;
    this->awgnGood = AWGN(0.0, 0.2, 1);
    this->awgnBad = AWGN(0.0, 1.5, 1);
}

Channel::~Channel() = default;


double Channel::get_random(double min, double max) {
    /* Returns a random double between min and max */
    return (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
}

std::vector<double> Channel::AWGNChannel(std::vector<int>& word){
    //std::cout << "word size: " << word.size() << std::endl;
    std::vector<double> received_word(word.size());
    std::vector<double> noiseVector = awgn.generateNoiseSamples();
    for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + noiseVector[i];
    }
    return received_word;
}  

std::pair<std::vector<double>, std::vector<double> > Channel::markovChannel(std::vector<double>& word){
    std::vector<double> received_word;
    // Vector that collects the different variance fo the symbols
    std::vector<double> differentVariance;
    // Every symbols of the modulated word has a different variance, chosen by the Markov channel
    // For the generation of the random number
    // Distribution in range [1, 10]
    // Initial state for the Markov chain

    for(int i = 0; i < word.size(); i++){
        double random_number = get_random(0.0, 1.0);
        std::vector<double> goodNoiseVector = awgnGood.generateNoiseSamples();
        std::vector<double> badNoiseVector = awgnBad.generateNoiseSamples();
        //std::cout << "random: " << random_number << std::endl;
        //std::cout << "state: " << actualState << std::endl;
        // The two states Markov chain
        if(random_number < 0.9999){
            if (actualState == 0){
                received_word.push_back(word[i] + goodNoiseVector[0]);
                differentVariance.push_back(0.2);
            }else{
                received_word.push_back(word[i] + badNoiseVector[0]);
                differentVariance.push_back(1.5);
            }
        }else{
            if (actualState == 0){
                received_word.push_back(word[i] + badNoiseVector[0]);
                differentVariance.push_back(1.5);
                actualState = 1;
            }else{
                received_word.push_back(word[i] + goodNoiseVector[0]);
                differentVariance.push_back(0.2);
                actualState = 0;
            } 
        }
    }
    return std::make_pair(received_word, differentVariance);
}

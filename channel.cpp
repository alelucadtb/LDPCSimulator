#include "channel.h"

Channel::Channel(AWGN& awgn) : awgn(awgn) {
}

Channel::Channel() {
    this->awgn = AWGN();
}

Channel::~Channel() = default;

std::vector<double> Channel::AWGNChannel(std::vector<int>& word){
    //std::cout << "word size: " << word.size() << std::endl;
    std::vector<double> received_word(word.size());
    std::vector<double> noiseVector = awgn.generateNoiseSamples();
    for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + noiseVector[i];
    }
    return received_word;
}  

std::pair<std::vector<double>, double> Channel::markovChannel(std::vector<double>& word){
    std::srand(std::time(0));
    std::vector<double> received_word(word.size());
    int random_number = rand() % 2;
    // Good state situation
    AWGN awgnGood(0.0, 0.1, word.size());
    // Bad state situation
    AWGN awgnBad(0.0, 2.0, word.size());
    std::vector<double> goodNoiseVector = awgnGood.generateNoiseSamples();
    std::vector<double> badNoiseVector = awgnBad.generateNoiseSamples();

    if(random_number == 0){
        for (int i = 0; i < word.size(); i++) {
            received_word[i] = word[i] + goodNoiseVector[i];
        }
        double variance = 0.1;
        return std::make_pair(received_word, variance);
    }else{
        for (int i = 0; i < word.size(); i++) {
            received_word[i] = word[i] + badNoiseVector[i];
        }
        double variance = 2.0;
        return std::make_pair(received_word, variance);
    }
}

std::tuple<std::vector<double>, double, std::vector<double>> Channel::storeMarkovChannel(std::vector<double>& word){
    std::srand(std::time(0));
    std::vector<double> received_word(word.size());
    int random_number = rand() % 2;
    // Good state situation
    AWGN awgnGood(0.0, 0.1, word.size());
    // Bad state situation
    AWGN awgnBad(0.0, 2.0, word.size());
    std::vector<double> goodNoiseVector = awgnGood.generateNoiseSamples();
    std::vector<double> badNoiseVector = awgnBad.generateNoiseSamples();

    if(random_number == 0){
        for (int i = 0; i < word.size(); i++) {
            received_word[i] = word[i] + goodNoiseVector[i];
        }
        double variance = 0.1;
        return std::make_tuple(received_word, variance, goodNoiseVector);
    }else{
        for (int i = 0; i < word.size(); i++) {
            received_word[i] = word[i] + badNoiseVector[i];
        }
        double variance = 2.0;
        return std::make_tuple(received_word, variance, badNoiseVector);
    }
}
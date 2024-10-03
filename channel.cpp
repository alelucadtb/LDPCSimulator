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
    for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + awgn.generateNoiseSamples()[i];
    }
    return received_word;
}  

std::vector<double> Channel::markovChannel(std::vector<double>& word){
    std::srand(std::time(0));
    std::vector<double> received_word(word.size());
    int random_number = rand() % 2;
    // Good state situation
    AWGN awgnGood(0.0, 0.25, word.size());
    // Bad state situation
    AWGN awgnBad(0.0, 1.0, word.size());

    if(random_number == 0){
        for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + awgnGood.generateNoiseSamples()[i];
        }
    }else{
        for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + awgnBad.generateNoiseSamples()[i];
        }
    }

    return received_word;
}
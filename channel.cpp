#include "channel.h"

Channel::Channel(AWGN& awgn) : awgn(awgn) {
}

Channel::~Channel() = default;

std::vector<double> Channel::AWGNChannel(std::vector<int>& word) {
    std::cout << "word size: " << word.size() << std::endl;
    std::vector<double> received_word(word.size());
    for (int i = 0; i < word.size(); i++) {
        received_word[i] = word[i] + awgn.generateNoiseSamples()[i];
    }
    return received_word;
}   
#ifndef CHANNEL_H
#define CHANNEL_H

#include "word.h"
#include "pam.h"
#include "AWGN.h"

class Channel {
    public:
    
    /*Constructor*/
    Channel(AWGN& awgn);
    Channel();

    /*Destructor*/
    ~Channel();

    /**
     * Representation of the AWGN channel
     * @param word: the word to be transmitted
     * @return the word after passing through the channel
    */
    std::vector<double> AWGNChannel(std::vector<int>& word);

    /**
     * Representation of the Markov model channel
     * @param word: the modulated word to be trasmitted
     * @return the word after passing through the channel
     * @return the different variances used by the channel
     */
    std::pair<std::vector<double>, std::vector<double> > markovChannel(std::vector<double> word);

    private:

    AWGN awgn;
    int actualState = 0;

};

#endif // CHANNEL_H
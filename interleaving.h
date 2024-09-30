#ifndef INTERLEAVING_H
#define INTERLEAVING_H

#include "word.h"
#include "pam.h"
#include "AWGN.h"

class Interleaving{
    public:

    /*Constructor*/
    Interleaving();

    /*Destructor*/
    ~Interleaving();

    std::vector<std::vector<double>> interleaving(std::vector<std::vector<int>> matrix);

    std::vector<std::vector<double>> deinterleaving(std::vector<std::vector<double>> matrix);


};


#endif // INTERLEAVING_H
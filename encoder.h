#ifndef ENCODER_H
#define ENCODER_H

#include <string>
#include "paritycheckmatrix.h"
#include "word.h"

class Encoder {
    public:

    /*Constructor*/ 
    Encoder(ParityCheckMatrix& matrix);

    /*Destructor*/
    ~Encoder();

    /**
     * Transforms the parity check matrix H into the generating matrix G
     * @return A vector of vectors representing the generating matrix G
    */
    std::vector<std::vector<int>> getGeneratingMatrix();

    /**
     * Encodes a word using the generating matrix G
     * @param word The word to encode
     * @return The encoded word
    */
    Word encode(Word& word);

    private:

    ParityCheckMatrix& matrix;

    /**
     * Applies Gaussian elimination to the parity check matrix H
     * @return A vector of vectors representing the generating matrix G
    */
    std::vector<std::vector<int>> gauss_reduce(Encoder& encoder);
};

#endif // ENCODER_H



#ifndef PAM_H
#define PAM_H

#include <vector>
#include "paritycheckmatrix.h"
#include "word.h"
#include "encoder.h"
#include "gray.h"
class PAM {
    public:

    /*Constructor*/
    PAM(Word& word, int M);

    /*Destructor*/
    ~PAM();

    private:
    /*The word to be modulated*/
    Word word;
    /*The number of symbols of the M-PAM*/
    int M;
    
    public:
    /*Get the number of symbols of the M-PAM*/
    int getM();
     /**
     * A method to do the PAM modulation
     * @param word The word to be modulated
     * @return The modulated word
    */    
    std::vector<int> PAMModulate(Word& word);

    /**
     * A method to do the M-PAM modulation
     * @param M The number of symbols of the M-PAM
     * @param word The word to be modulated
     * @return The modulated word
    */
    std::vector<int> MPAMModulate(Word& word);
};

#endif // PAM_H 
#ifndef ERROR_H
#define ERROR_H

#include <vector>
#include <cmath>
#include <iostream>
#include "word.h"

class Error{
    public:
        /*Constructor and Destructor*/
        Error();
        ~Error();

        double calculateError(std::vector<int> &received, Word &original);
};

#endif
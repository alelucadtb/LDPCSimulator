/*
    This file contains the functions to generate the parity check matrix following the IEEE 802.11-2020 standard
*/

#ifndef PARITYCHECKMATRIX_H
#define PARITYCHECKMATRIX_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "word.h"

using namespace std;

class ParityCheckMatrix {
    public:

    int N;
    int K;

    /*Constructor*/     
    ParityCheckMatrix(int n, int z, double R, const std::string& filename);
    /*Destructor*/
    ~ParityCheckMatrix();

    /**
     * Returns the base matrix
     * @return A vector of vectors representing the base matrix
    */
    std::vector<std::vector<int> > getBaseMatrix() const {
        return base_matrix;
    }

    /**
     * Generates a parity check matrix in a binary form
     * @return A vector of vectors representing the generated matrix
    */
    std::vector<std::vector<int> > getBinaryMatrix() const;

    /**
     * Returns the number of rows N-K
     * @return The number of rows N-K
    */
    int getRows() const {
        return (n / z) - (R * (n / z));
    }

    /**
     * Returns the number of columns N
     * @return The number of columns N
    */
    int getCols() const {
        return n / z;
    }

    int getBitCols() const {
        return n;
    }

    int getBitRows() const {
        return n - (R * n);
    }

    /**
     * Loads a base matrix from a file
     * @param filename The name of the file to read from
     * @return A vector of vectors representing the base matrix
    */
    std::vector<std::vector<int> > loadBaseMatrix(const std::string& filename);

    /**
     * Writes the binary matrix to a file
     * @param filename The name of the file to write to
    */
    void writeBinaryMatrix(const std::string& filename) const;

    /**
         * Verifies if a vector is a codeword
         * @param word The word to verify
         * @return True if the word is a codeword, false otherwise
    */
    bool isCodeword(Word& word) const;

    /**
         * Verifies if a vector is a codeword
         * @param word The word to verify
         * @return True if the word is a codeword, false otherwise
    */
    bool isCodewordVector(std::vector<int>& word) const;
    
    private:

    int n;
    int z;
    double R;
    std::vector<std::vector<int> > base_matrix;

    /**
     * Creates an identity matrix P_i and performs a cyclic permutation of i positions
     * @param size The size of the square matrix
     * @param i The number of positions to shift
     * @return A vector of vectors representing the permuted identity matrix
    */
    std::vector<std::vector<int> > create_permuted_identity_matrix(std::size_t size, std::size_t i) const;

    /**
     * Creates a null matrix of size z x z
     * @param z The size of the square matrix
     * @return A vector of vectors representing the null matrix
    */
    std::vector<std::vector<int> > create_null_matrix(std::size_t z) const;

};

#endif // PARITYCHECKMATRIX_H
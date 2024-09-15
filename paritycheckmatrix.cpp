#include "paritycheckmatrix.h"
#include "word.h"
#include <vector>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

/*Constructor*/
ParityCheckMatrix::ParityCheckMatrix(int n, int z, double R, const std::string& filename) { 
    this->base_matrix = loadBaseMatrix(filename);
    this->n = n;
    this->z = z;
    this->R = R;
}

/*Destructor*/
ParityCheckMatrix::~ParityCheckMatrix() {
    // Destructor
}   

std::vector<std::vector<int>> ParityCheckMatrix::create_permuted_identity_matrix(std::size_t size, std::size_t i) const {
    std::vector<std::vector<int>> matrix(size, std::vector<int>(size, 0));
    for (std::size_t row = 0; row < size; ++row) {
        // Perform a cyclic permutation of i positions
        matrix[row][(row + i) % size] = 1;  
    }
    return matrix;
}

std::vector<std::vector<int>> ParityCheckMatrix::create_null_matrix(std::size_t z) const {
    return std::vector<std::vector<int>>(z, std::vector<int>(z, 0));
}

std::vector<std::vector<int>> ParityCheckMatrix::getBinaryMatrix() const {
    int N = n / z;
    int K = R * N;

    // Check if the base matrix has the correct dimensions
    if (base_matrix.size() != (N - K) || base_matrix[0].size() != N) {
        throw std::invalid_argument("The dimension of the input matrix must be (N-K)xN");
    }

    // Initialize the generated matrix with the correct dimensions  
    std::vector<std::vector<int>> generated_matrix((N - K) * z, std::vector<int>(N * z, 0));

    // Fill the generated matrix with the appropriate sub-matrices
    for (std::size_t row = 0; row < base_matrix.size(); ++row) {
        for (std::size_t col = 0; col < base_matrix[row].size(); ++col) {
            std::vector<std::vector<int>> sub_matrix;
            if (base_matrix[row][col] == -1) {
                sub_matrix = create_null_matrix(z);
            } 
            else {
                sub_matrix = create_permuted_identity_matrix(z, base_matrix[row][col]);
            }

            // Copy the sub-matrix into the generated matrix
            for (std::size_t sub_row = 0; sub_row < z; ++sub_row) {
                for (std::size_t sub_col = 0; sub_col < z; ++sub_col) {
                    generated_matrix[row * z + sub_row][col * z + sub_col] = sub_matrix[sub_row][sub_col];
                }
            }
        }
    }
    return generated_matrix;
}

std::vector<std::vector<int>> ParityCheckMatrix::loadBaseMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<int>> matrix;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::vector<int> row;
        int value;
        while (ss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

void ParityCheckMatrix::writeBinaryMatrix(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<int>> binary_matrix = getBinaryMatrix();
    for (const auto& row : binary_matrix) {
        for (std::size_t col = 0; col < row.size(); ++col) {
            file << row[col];
            if (col < row.size() - 1) {
                file << "\t";
            }
        }
        file << "\n";
    }

    file.close();
}

bool ParityCheckMatrix::isCodeword(Word& word) const{
    std::vector<std::vector<int>> binary_matrix = getBinaryMatrix();
    std::vector<int> result(binary_matrix.size(), 0);
    for (std::size_t i = 0; i < binary_matrix.size(); ++i) {
        for (std::size_t j = 0; j < binary_matrix[i].size(); ++j) {
            result[i] ^= (binary_matrix[i][j] * word.get(j));
        }
    }
    for (const auto& val : result) {
        if (val != 0) {
            return false;
        }
    }
    return true;
}

bool ParityCheckMatrix::isCodewordVector(std::vector<int>& word) const{
    std::vector<std::vector<int>> binary_matrix = getBinaryMatrix();
    std::vector<int> result(binary_matrix.size(), 0);
    for (std::size_t i = 0; i < binary_matrix.size(); ++i) {
        for (std::size_t j = 0; j < binary_matrix[i].size(); ++j) {
            result[i] ^= (binary_matrix[i][j] * word[j]);
        }
    }
    for (const auto& val : result) {
        if (val != 0) {
            return false;
        }
    }
    return true;
}
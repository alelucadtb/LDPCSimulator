#include "encoder.h"

/*Constructor*/
Encoder::Encoder(ParityCheckMatrix& matrix) : matrix(matrix) {
}

/*Destructor*/
Encoder::~Encoder() = default;

std::vector<std::vector<int>> Encoder::gauss_reduce(Encoder& encoder) {
    int rows = matrix.getBitRows();
    int cols = matrix.getBitCols();
    std::vector<std::vector<int>> H = matrix.getBinaryMatrix();

    for (int i = 0; i < rows; ++i) {
        // Make the diagonal contain all 1's
        if (H[i][i + (cols - rows)] == 0) {
            for (int j = i + 1; j < rows; ++j) {
                if (H[j][i + (cols - rows)] == 1) {
                    std::swap(H[i], H[j]);
                    break;
                }
            }
        }

        // Make all rows below this one 0 in current column
        for (int j = i + 1; j < rows; ++j) {
            if (H[j][i + (cols - rows)] == 1) {
                for (int k = 0; k < cols; ++k) {
                    H[j][k] ^= H[i][k];
                }
            }
        }
    }

    // Make all rows above this one 0 in current column
    for (int i = rows - 1; i >= 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            if (H[j][i + (cols - rows)] == 1) {
                for (int k = 0; k < cols; ++k) {
                    H[j][k] ^= H[i][k];
                }
            }
        }
    }

    return H;
}

std::vector<std::vector<int>> Encoder::getGeneratingMatrix() {
    std::vector<std::vector<int>> H = gauss_reduce(*this);
    int n = matrix.getBitCols();
    int k = n - matrix.getBitRows();
    std::vector<std::vector<int>> G(n, std::vector<int>(k, 0));

    // Fill the first k x k elements of G with an identity matrix
    for (int i = 0; i < k; ++i) {
        G[i][i] = 1;
    }

    // Fill the last (n-k) x k elements of G with the first (n-k) x k elements of H
    for (int i = 0; i < n - k; ++i) {
        for (int j = 0; j < k; ++j) {
            G[k + i][j] = H[i][j];
        }
    }
    return G;
    
}

Word Encoder::encode(Word& word) {
    std::vector<std::vector<int>> G = getGeneratingMatrix();
    int k = G[0].size();
    int n = G.size();
    std::vector<int> encoded_word(n, 0);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            encoded_word[i] ^= (G[i][j] * word.get(j));
        }
    }

    return Word(encoded_word);
}
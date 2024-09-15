#include "pam.h"
#include <cmath>

/*Constructor*/
PAM::PAM(Word& word, int M) : word(word), M(M) {
}

/*Destructor*/
PAM::~PAM() = default;

int PAM::getM(){
    return M;
}

std::vector<int> PAM::PAMModulate(Word& word) {
    std::vector<int> modulated_data(word.length());
    for (int i = 0; i < word.length(); ++i) {
        if (word.get(i) == 1) {
            modulated_data[i] = 1;
        } else {
            modulated_data[i] = -1;
        }
    }
    return modulated_data;
}

std::vector<int> PAM::MPAMModulate(Word& word) {
    std::vector<int> modulated_data;
    Gray gray;
    /*Get the number of bits per symbol*/
    int n = log2(M);
    std::string gray_total = "";
    /*Generate the gray codes*/
    std::vector<std::string> gray_codes = gray.generateGrayarr(n);
    /*Generate the gray codes in integer*/
    std::vector<int> gray_codes_int(gray_codes.size());
    for (int i = 0; i < gray_codes.size(); i++) {
       gray_codes_int[i] = 2*(i)-(M-1);
    }
    /*Modulate the data*/
    for (int i = 0; i < word.length(); i = i + n) {
        gray_total = "";
        /*Convert the bits to a gray code*/
        for (int j = 0; j < n; j++) {
            int value = word.get(i + j);
            std::string value_str = std::to_string(value);
            gray_total = gray_total + value_str;
        }
        /*Find the corresponding gray code*/
        for (int j = 0; j < gray_codes.size(); j++) {
            if (gray_codes[j] == gray_total) {
                modulated_data.push_back(gray_codes_int[j]);
                break;
            }
        }
    }
    return modulated_data;
}
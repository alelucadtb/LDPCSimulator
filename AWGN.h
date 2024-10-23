#ifndef AWGN_h
#define AWGN_h

#include <iostream>
#include <random>
#include <vector>

class AWGN {
    
public:
    AWGN(); // Default constructor
    AWGN(double mean, double variance, int numberOfSamples); // Constructor with initializers
    AWGN(double SNR, int numberOfSamples); // Constructor with initializers (SNR)
    
    std::vector<double> generateNoiseSamples();
    
private:
    double mean;
    double variance;
    double stddev;
    int numberOfSamples;
    double sigma;
    bool isSNRMode;
    
};

#endif /* AWGN_h */
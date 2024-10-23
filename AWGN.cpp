#include "AWGN.h"

AWGN::AWGN (){
    // default values are assigned
    this->mean = 0.0;
    this->variance = 1.0;
    this->numberOfSamples = 100;
    this->isSNRMode = false;  
    //std::cerr << "WARNING!, using default initializer is potentially dangerous, default parameters values are used!" << std::endl;
}

AWGN::AWGN (double mean, double variance, int numberOfSamples){
    this->mean = mean;
    this->variance = variance;
    this->stddev = sqrt(variance);
    this->numberOfSamples = numberOfSamples;
    this->isSNRMode = false;
}

AWGN::AWGN (double SNR, int numberOfSamples){
    this->mean = 0.0;
    this->variance = 1.0;
    this->stddev = sqrt(variance);
    this->numberOfSamples = numberOfSamples;
    this->sigma = sqrt(pow(10,(-SNR/10)));
    this->isSNRMode = true;
    std::cout << "WARNING! singal power is assumed to be normalized" << std::endl;  
}

std::vector<double> AWGN::generateNoiseSamples(){
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normalDistribution(this->mean,this->stddev);
    
    // Dynamically allocate noise samples
    std::vector<double> samples(this->numberOfSamples);
    
    if (this->isSNRMode) {
        // calculate with mean and variance
        for (int i = 0; i<this->numberOfSamples; i++) {
            samples[i] = (this->sigma) * normalDistribution(gen);
        }
        
    } else {
        // calculate with SNR
        for (int i = 0; i<this->numberOfSamples; i++) {
            samples[i] = normalDistribution(gen);
        }
        
    }
    
    return samples;
}
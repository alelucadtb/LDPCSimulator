// C++ program to generate n-bit Gray codes
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Gray {
public:
    /*Constructor*/
    Gray();
    /*This function generates all n bit Gray codes and prints the generated codes*/                    
    std::vector<std::string> generateGrayarr(int n);
};  
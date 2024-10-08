#include "gray.h"
using namespace std;

// C++ program to generate n-bit Gray codes
#include <iostream>
#include <string>
#include <vector>
using namespace std;

Gray::Gray(){

}
// This function generates all n bit Gray codes and prints the
// generated codes
std::vector<std::string> Gray::generateGrayarr(int n)
{
    // base case
    if (n <= 0)
        return std::vector<std::string>();
 
    // 'arr' will store all generated codes
    std::vector<std::string> arr;
 
    // start with one-bit pattern
    arr.push_back("0");
    arr.push_back("1");
 
    // Every iteration of this loop generates 2*i codes from previously
    // generated i codes.
    int i, j;
    for (i = 2; i < (1<<n); i = i<<1)
    {
        // Enter the previously generated codes again in arr[] in reverse
        // order. Nor arr[] has double number of codes.
        for (j = i-1 ; j >= 0 ; j--)
            arr.push_back(arr[j]);
 
        // append 0 to the first half
        for (j = 0 ; j < i ; j++)
            arr[j] = "0" + arr[j];
 
        // append 1 to the second half
        for (j = i ; j < 2*i ; j++)
            arr[j] = "1" + arr[j];
    }
    return arr;
}
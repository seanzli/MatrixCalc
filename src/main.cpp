#include "../include/matrix.hpp"

#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
    // 
    double test[4] = {1,2,3,4};
    MatrixCalc test_m(test, 2, 2);

    
    std::cout << MatrixCalc::eye(3, 5) << std::endl;
    

    return 0;
}


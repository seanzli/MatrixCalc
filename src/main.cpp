#include "../include/matrix.hpp"

#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
    // 
    double test[4] = {1,2,3,4};
    MatrixCalc test_m(test, 2, 2);

    std::cout << "orj :" << test_m << std::endl;
    std::cout << "sum :" << test_m.sum() << std::endl;
    std::cout << "sum_row :" << test_m.sum_row() << std::endl;
    std::cout << "sum_col :" << test_m.sum_col() << std::endl;
    std::cout << "mean :" << test_m.mean() << std::endl;
    std::cout << "mean_row :" << test_m.mean_row() << std::endl;
    std::cout << "mean_col :" << test_m.mean_col() << std::endl;
    

    return 0;
}



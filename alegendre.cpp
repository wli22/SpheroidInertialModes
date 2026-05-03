#include <iostream>
#include <cmath>
#include <iomanip>
#include "ALegendre.hpp"


int main() {

    int l = 1000;
    int m = 500;
    double theta = 0.5;
    double x = std::cos(theta);

    ALegendre<double> alegendre;

    double plm = alegendre.ALegendre_compute(l, m, x);
    
    std::cout << "l = " << l << "   m = " << m << "   plm = " << std::scientific << std::setprecision(15) << plm << "\n";

    return 0;
}
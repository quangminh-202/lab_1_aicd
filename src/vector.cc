#include <vector/vector.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <complex>
using namespace std;


int main() {
    // Example usage for a vector of doubles
    double double_coords[] = { 3.0, 4.0 };
    my_vector::vector<double> double_vector(2, double_coords);
    my_vector::vector<double> result_double = calculate_perpendicular_unit_vector(double_vector);
    std::cout << "Perpendicular unit vector for double vector: " << result_double << std::endl;

    // Example usage for a vector of complex numbers
   /* std::complex<double> complex_coords[] = { std::complex<double>(1.0, 2.0), std::complex<double>(3.0, 4.0) };
    my_vector::vector<std::complex<double>> complex_vector(2, complex_coords);
    my_vector::vector<std::complex<double>> result_complex = calculate_perpendicular_unit_vector(complex_vector);
    std::cout << "Perpendicular unit vector for complex vector: " << result_complex << std::endl;*/

    return 0;
}

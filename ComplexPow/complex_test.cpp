#include <iostream>
#include <string>
#include <math.h>
#include "imgry.h"

template <class T>
void PrintComplex(Imaggi<T>& z){
    T x_real = z.GetRealPart();
    T y_imag = z.GetImagPart();

    std::cout << x_real <<" + " << y_imag << "i" << std::endl;
}

template <class T>
void PrintComplex(Imaggi<T>& z, const std::string& w){
    T x_real = z.GetRealPart();
    T y_imag = z.GetImagPart();

    std::cout << w << " = " << x_real <<" + " << y_imag << "i" << std::endl;
}

int main (){

    double x1 = 3.4;
    double y1 = 5;
    Imaggi<double> z1(x1,y1);
    double x2 = 5.1;
    double y2 = 2;
    Imaggi<double> z2(x2,y2);

    PrintComplex(z1, "z1");
    PrintComplex(z2, "z2");

    std::cout << "Conjugado z1" << std::endl;
    Imaggi<double> conj = Conjugate(z1);
    PrintComplex(conj, "z1*");

    std::cout << "Sum z1 + z2" << std::endl;
    Imaggi<double> sum = z1 + z2;
    PrintComplex(sum, "z1+z2");

    std::cout << "Product z1 * z2" << std::endl;
    Imaggi<double> product = z1 * z2;
    PrintComplex(product, "z1*z2");

    /*std::cout << "Division z1 / z2" << std::endl;
    Imaggi<double> div = z1 / z2;
    PrintComplex(div, "z1/z2");*/

    std::cout << "Exponentiation z2 ^ 3" << std::endl;
    Imaggi<double> power = z2 ^ 3;
    PrintComplex(power, "z1 ^ 3");

    return 0;
}

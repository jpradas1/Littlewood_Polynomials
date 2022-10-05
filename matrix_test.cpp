#include <iostream>
#include <cmath>
#include "matrix.h"

int main(){

    // double a = 2.4;
    // double Data1[9] = {1, 0, 1, 2, -1, 0, 3, 0, 2};
    // double Data2[12] = {0, 1, 5, 3, 2, 4, 2, -9, 3, 1, 6, 2};
    // double Data3[9] = {8, 0, 5, 2, 11, 4, 3, 7, 2};
    // nmMatrix data1(3,3,Data1);
    // nmMatrix data2(3,4,Data2);
    // nmMatrix data3(3,3,Data3);
    // nmMatrix identity(20,20,'i');
    nmMatrix data4(5,5);
    data4.uniform_random();
    data4.showMatrix();
    nmMatrix data5(5,5);
    data5.uniform_random();
    data5.showMatrix();

    /* std::cout << "original matrix" << std::endl;
    PrintMatrix(data1);
    std::cout << "\n";
    PrintMatrix(data2);
    std::cout << "\n";
    PrintMatrix(data3);
    std::cout << "\n";
    PrintMatrix(data4); */

    /* std::cout << "suma data1+data3" << std::endl;
    nmMatrix sum13 = data1 - data3;
    PrintMatrix(sum13);
    std::cout << "suma data1+data2" << std::endl;
    nmMatrix sum23 = data1 - data2;
    PrintMatrix(sum23);
    std::cout << "suma data1+=data3" << std::endl;
    data1 -= data3;
    PrintMatrix(data1);
    std::cout << "\n";
    PrintMatrix(data3); */

   /*  std::cout << "original matrix" << std::endl;
    PrintMatrix(data1);
    std::cout << "\n";
    PrintMatrix(data2);
    std::cout << "\n";
    PrintMatrix(data3); */


    /* std::cout << "multiplication" << std::endl;
    nmMatrix scalarL = a * data3;
    PrintMatrix(scalarL);
    std::cout << "\n";
    nmMatrix scalarR = data3 * a;
    PrintMatrix(scalarR);
    std::cout << "\n";
    nmMatrix product12 = data1 * data2;
    PrintMatrix(product12); */

   /*  if(data1 == data2)
        std::cout << "Peguelo perro" << std::endl;
    else if(data1 == data4)
        std::cout << "Buena la rata" << std::endl; */

    return 0;
}

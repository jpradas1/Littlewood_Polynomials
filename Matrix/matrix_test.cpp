#include <iostream>
#include <string>
#include <math.h>
#include "matrix.h"

template<class T>
void PrintMatrix(nmMatrixN<T>& nmMatrix){
    int rows = nmMatrix.GetNumRows();
    int cols = nmMatrix.GetNumCols();

    for(int ii = 0 ; ii<  rows; ii++){
        for(int jj = 0 ; jj < cols; jj++){
            std::cout << nmMatrix.GetElements(ii, jj) << "  ";
        }
        std::cout << "\n";
    }
}

int main(){

    int Data1[9] = {1, 0, 1, 2, -1, 0, 3, 0, 2};
    int Data2[12] = {0, 1, 5, 3, 2, 4, 2, -9, 3, 1, 6, 2};
    int Data3[9] = {8, 0, 5, 2, 11, 4, 3, 7, 2};
    nmMatrixN<int> data1(3,3,Data1);
    nmMatrixN<int> data2(3,4,Data2);
    nmMatrixN<int> data3(3,3,Data3);

    std::cout << "original matrix" << std::endl;
    PrintMatrix(data1);
    std::cout << "\n";
    PrintMatrix(data2);
    std::cout << "\n";
    PrintMatrix(data3);

    std::cout << "suma data1+data3" << std::endl;
    nmMatrixN<int> sum13 = data1 + data3;
    PrintMatrix(sum13);

    std::cout << "multiplication data1*data2" << std::endl;
    nmMatrixN<int> product12 = data1 * data2;
    PrintMatrix(product12);

    return 0;
}

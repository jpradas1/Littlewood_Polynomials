#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>

template <class T>
class nmMatrixN
{
    public:
        nmMatrixN ();//Define an array with only one elemente == 0
        //Define an array with n rows ans m cols where its elements are zeros
        nmMatrixN (int nRows, int mCols);
        //Let an input array, we fill a new T object with the array's original elements
        nmMatrixN (int nRows, int mCols, const T *arrayData);
        //Make a copy of the nmMatrix
        nmMatrixN (const nmMatrixN<T>& nmMatrix);

        //And obviously the destructor
        ~nmMatrixN ();

        //Wipe out the private matrixData object and resize it
        bool resize (int newRows, int newCols);
        //Obtain elements in the nth row and mth column
        T GetElements (int row, int col);
        //Here we can define each matrix's component element by element
        bool EbyElement (int row, int col, T element);

        int GetNumRows(); int GetNumCols();// Name literally says what they do

        //it's time to determinate basic operation among matrices (+,-,*,=)


    private:
        //scanner: let nth row and mth columns of a matrix gives me position into array
        int GetPosition (int rows, int cols);
    private:
        T *matrixData; //object where it saves information
        int N_Rows, M_Cols, NM_Elements;

};
// #######################   Constructors   #######################

//Default matrix
template <class T>
nmMatrixN<T>::nmMatrixN (){
    N_Rows = 1;
    M_Cols = 1;
    NM_Elements = 1;
    matrixData = new T [NM_Elements];
    matrixData[0] = 0.0;
}

//Matrix(n*m) initialized with just zeros
template <class T>
nmMatrixN<T>::nmMatrixN (int nRows, int mCols){
    N_Rows = nRows;
    M_Cols = mCols;
    NM_Elements = N_Rows*M_Cols;
    matrixData = new T [NM_Elements];
    for (unsigned int ii = 0; ii < NM_Elements ; ii++){
        matrixData[ii] = 0.0;
    }
}

//Matrix(n*m) filled with respect array's elements
template <class T>
nmMatrixN<T>::nmMatrixN (int nRows, int mCols, const T *arrayData){
    N_Rows = nRows;
    M_Cols = mCols;
    NM_Elements = N_Rows*M_Cols;
    matrixData = new T [NM_Elements];
    for (unsigned int ii = 0; ii < NM_Elements ; ii++){
        matrixData[ii] = arrayData[ii];
    }
}

//Copy matrix nmMatrixN like
template <class T>
nmMatrixN<T>::nmMatrixN (const nmMatrixN<T>& nmMatrix){
    N_Rows = nmMatrix.N_Rows;
    M_Cols = nmMatrix.M_Cols;
    NM_Elements = nmMatrix.NM_Elements;
    matrixData = new T [NM_Elements];
    for (unsigned int ii = 0; ii < NM_Elements ; ii++){
        matrixData[ii] = nmMatrix.matrixData[ii];
    }
}

//Destructor doesn't allow memory leak
nmMatrixN<T>::~nmMatrixN (){
    if(matrixData != nullptr){
        delete [] matrixData;
    }
}

// ####################### Some elemental funcrions #######################

template <class T>
bool resize (int newRows, int newCols){
    N_Rows = newRows;
    M_Cols = newCols;
    NM_Elements = N_Rows*M_Cols;
    delete[] matrixData;
    matrixData = new T [NM_Elements]
    if (matrixData != nullptr){
        for(unsigned int ii = 0 ; ii < NM_Elements ; ii++){
            da
        }
    }
}

#endif

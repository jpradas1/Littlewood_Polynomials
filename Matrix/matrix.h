#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>

template <class T>
class nmMatrixN
{
    public:
        nmMatrixN ();//Define an array with only one elemente == 0
        //Define an array with n rows and m cols where its elements are zeros
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

        //it's time to determinate basic operation on matrices (+,-,*,=)
        //by overloading these operators

        bool operator== (const nmMatrixN<T>& nmMatrix);

        //The following functions are defined through 'friend' pointer which helps
        //to define functions that don't belong this class but need the private parameters
        //and perhaps another instances

        template <class S> friend nmMatrixN<S> operator+ (const nmMatrixN<S>& lMatrix, const nmMatrixN<S>& rMatrix); //sum on two matrix
        template <class S> friend nmMatrixN<S> operator+ (const S& lscaler, const nmMatrixN<S>& rMatrix); //scaler + matrix
        template <class S> friend nmMatrixN<S> operator+ (const nmMatrixN<S>& lMatrix, const S& rscaler); //matrix + scaler

        template <class S> friend nmMatrixN<S> operator- (const nmMatrixN<S>& lMatrix, const nmMatrixN<S>& rMatrix); //subtraction on two matrix
        template <class S> friend nmMatrixN<S> operator- (const S& lscaler, const nmMatrixN<S>& rMatrix);
        template <class S> friend nmMatrixN<S> operator- (const nmMatrixN<S>& lMatrix, const S& rscaler);

        template <class S> friend nmMatrixN<S> operator* (const nmMatrixN<S>& lMatrix, const nmMatrixN<S>& rMatrix); //product on two matrix
        template <class S> friend nmMatrixN<S> operator* (const S& lscaler, const nmMatrixN<S>& rMatrix); //product by scaler lefthandside
        template <class S> friend nmMatrixN<S> operator* (const nmMatrixN<S>& lMatrix, const S& rscaler); //product by scaler righthandside


    private:
        //scanner: let nth row and mth column of a matrix gives me position into array
        int GetPosition (int row, int col);
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
    for (int ii = 0; ii < NM_Elements ; ii++){
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
    for (int ii = 0; ii < NM_Elements ; ii++){
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
    for (int ii = 0; ii < NM_Elements ; ii++){
        matrixData[ii] = nmMatrix.matrixData[ii];
    }
}

//Destructor doesn't allow memory leak
nmMatrixN<T>::~nmMatrixN (){
    if(matrixData != nullptr){
        delete [] matrixData;
    }
}

// ####################### Some elemental functions #######################

//wipe out and generate a new matrix with different dimension
template <class T>
bool nmMatrixN<T>::resize (int newRows, int newCols){
    N_Rows = newRows;
    M_Cols = newCols;
    NM_Elements = N_Rows*M_Cols;
    delete[] matrixData;
    matrixData = new T [NM_Elements];
    if (matrixData != nullptr){
        for(int ii = 0 ; ii < NM_Elements ; ii++){
            matrixData[ii] = 0.0; }

        return true;
    }
    else { return false;}
}

//Let certain index positions in matrix we obtain elements on that position
template <class T>
T nmMatrixN<T>::GetElements (int row, int col){
    int position_e = GetPosition(row, col);
    if ( position_e >= 0){
        return matrixData[position_e];
    }
    else return 0.0;
}

//Let certain index position on matrix, we obtain its respect position into array
template <class T>
int nmMatrixN<T>::GetPosition (int row, int col){
    if ((row >= 0) && (row < N_Rows) && (col >= 0) && (col < M_Cols)){
        return row*M_Cols + col;
    }
    else return -1;
}

//This allow to define a matrix element by element or replace a certain element by other
template <class T>
bool nmMatrixN<T>::EbyElement (int row, int col, T value){
    int position_e = GetPosition (row, col);
    if ( position_e >= 0 ){
        matrixData [position_e] = value;
        return true;
    }
    else return false;
}

template <class T>
int nmMatrixN<T>::GetNumRows(){
    return N_Rows;
}

template <class T>
int nmMatrixN<T>::GetNumCols(){
    return M_Cols;
}

// ####################### Operations on Matrix #######################

// SUM
template <class T>
nmMatrixN<T> operator+ (const nmMatrixN<T>& lMatrix, const nmMatrixN<T>& rMatrix){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    T *partial_res = new T[lElements];
    //int rRows = rMatrix.N_Rows;
    //int rCols = rMatrix.M_Cols;
    //int rElements = lRows*lCols;
    //if((lRows == rRows) && (lCols == rCols))
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix[ii] + rMatrix[ii];}

    nmMatrixN<T> result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
    //else return nmMatrixN<T> result(1,1);
}

template <class T>
nmMatrixN<T> operator+ (const T& lscaler, const nmMatrixN<T>& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows*rCols;
    T *partial_res = new T[rElements];
    for(int ii = 0; ii < rElements; ii++){
        partial_res[ii] = lscaler + rMatrix[ii];}

    nmMatrixN<T> result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}

template <class T>
nmMatrixN<T> operator+ (const nmMatrixN<T>& lMatrix, const T& rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    T *partial_res = new T[lElements];
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix[ii] + rscaler;}

    nmMatrixN<T> result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

// SUBTRACTION (which is essentially the same)

template <class T>
nmMatrixN<T> operator- (const nmMatrixN<T>& lMatrix, const nmMatrixN<T>& rMatrix){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    T *partial_res = new T[lElements];
    //int rRows = rMatrix.N_Rows;
    //int rCols = rMatrix.M_Cols;
    //int rElements = lRows*lCols;
    //if((lRows == rRows) && (lCols == rCols))
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix[ii] - rMatrix[ii];}
    nmMatrixN<T> result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
    //else return nmMatrixN<T> result(1,1);
}

template <class T>
nmMatrixN<T> operator- (const T& lscaler, const nmMatrixN<T>& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows*rCols;
    T *partial_res = new T[rElements];
    for(int ii = 0; ii < rElements; ii++){
        partial_res[ii] = lscaler - rMatrix[ii];}

    nmMatrixN<T> result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}

template <class T>
nmMatrixN<T> operator- (const nmMatrixN<T>& lMatrix, const T& rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    T *partial_res = new T[lElements];
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix[ii] - rscaler;}

    nmMatrixN<T> result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

//MULTIPLICATION (maybe longer forward I could implement multiplication by blocking method)

template <class T>
nmMatrixN<T> operator* (const nmMatrixN<T>& lMatrix, const nmMatrixN<T>& rMatrix){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;

    if (lCols == rRows){
        T *partial_res = new T[lRows*rCols];
        for(int ii = 0; ii < lRows; ii++){
            for(int jj = 0 ; jj < rCols; jj++){
                partial_res[ii*rCols+jj] = 0.0;
                for(int kk = 0; kk < lCols; kk++){
                 partial_res[ii*rCols+jj] += lMatrix[ii*lCols+kk] * rMatrix[kk*rCols+ii];
                }
            }
        }

        nmMatrixN<T> result(lRows, rCols, partial_res);
        delete[] partial_res;
        return result;
    }
    else {nmMatrixN<T> result(1,1);
        return result;}
}

template <class T>
nmMatrixN<T> operator* (const T& lscaler, const nmMatrixN<T>& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows * rCols;
    T *partial_res = new T[rElements];
    for(int ii = 0 ; ii < rElements; ii++){
        partial_res[ii] = lscaler * rMatrix[ii];
    }

    nmMatrixN<T> result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}

template <class T>
nmMatrixN<T> operator* (const nmMatrixN<T>& lMatrix, const T& rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows * lCols;
    T *partial_res = new T[lElements];
    for(int ii = 0 ; ii < lElements; ii++){
        partial_res[ii] = lMatrix[ii] * rscaler;
    }

    nmMatrixN<T> result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

//EQUALITY

template <class T>
bool nmMatrixN<T>::operator== (const nmMatrixN<T>& nmMatrix){
    //Checking coincidence matrices dimension
    if((this->N_Rows != nmMatrix.N_Rows) && (this->M_Cols != nmMatrix.M_Cols)){
        return false;}
    //Checking coincidence element by element
    bool bug = true;
    for(int ii = 0; ii < this->NM_Elements ; ii++){
        if(this->matrixData[ii] != nmMatrix.matrixData[ii]){
            bug = false;
        }
    }

    return bug;
}

#endif

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include "./random/Random64.h"

Crandom rand64(1);

class nmMatrix
{
    public:
        nmMatrix ();//Define an array with only one elemente = 0
        //Define an array with n rows and m columns where its elements are zeros
        nmMatrix (int nRows, int mCols);
        //Let an input array, we fill a matrix using its data
        nmMatrix (int nRows, int mCols, const double * arrayData);
        //Make a copy of the nmMatrix object
        nmMatrix (const nmMatrix & nm_Matrix);
        //Construct the identity matrix
        nmMatrix (int nRows, int mCols, char identity);

        ~nmMatrix ();

        void resize (int newRows, int newCols);
        double GetElements (int row, int col);
        void replace (int row, int col, const double value);
        int GetNumRows(){return N_Rows;};
        int GetNumCols(){return M_Cols;};
        void showMatrix();

        // OPERATIONS

        bool operator== (const nmMatrix& nmMatrix);

        friend nmMatrix operator+ (const nmMatrix& lMatrix, const nmMatrix& rMatrix);
        friend nmMatrix operator+ (const double & lscaler, const nmMatrix& rMatrix);
        friend nmMatrix operator+ (const nmMatrix& lMatrix, const double & rscaler);

        friend nmMatrix operator+= (const nmMatrix& lMatrix, const nmMatrix& rMatrix);
        friend nmMatrix operator+= (const nmMatrix& lMatrix, const double & rscaler);

        friend nmMatrix operator- (const nmMatrix& lMatrix, const nmMatrix& rMatrix);
        friend nmMatrix operator- (const double & lscaler, const nmMatrix& rMatrix);
        friend nmMatrix operator- (const nmMatrix& lMatrix, const double & rscaler);

        friend nmMatrix operator-= (const nmMatrix& lMatrix, const nmMatrix& rMatrix);
        friend nmMatrix operator-= (const nmMatrix& lMatrix, const double & rscaler);

        friend nmMatrix operator* (const nmMatrix& lMatrix, const nmMatrix& rMatrix);
        friend nmMatrix operator* (const double & lscaler, const nmMatrix& rMatrix);
        friend nmMatrix operator* (const nmMatrix& lMatrix, const double & rscaler);

    public: //special functions
        //fill the matrix with random values between (0,1)
        void uniform_random();
        //fill the matrix using a normal distribution
        void normal_random(const double mu, const double sigma);

    private:
        //scanner: let nth row and mth column of a matrix gives the position into array
        int GetPosition (int row, int col);
    private:
        double * matrixData;
        int N_Rows, M_Cols, NM_Elements;
};
// %%%%%%%%%%%%%%%%%%%%%%%%%  Constructors  %%%%%%%%%%%%%%%%%%%%%%%%%


nmMatrix::nmMatrix (){
    N_Rows = 1;
    M_Cols = 1;
    NM_Elements = N_Rows * M_Cols;
    matrixData = new double [NM_Elements];
    matrixData[0] = 0.0;
}


nmMatrix::nmMatrix (int nRows, int mCols){
    N_Rows = nRows;
    M_Cols = mCols;
    NM_Elements = N_Rows * M_Cols;
    matrixData = new double [NM_Elements];
    for (int ii = 0; ii < NM_Elements ; ii++){
        *(matrixData + ii) = 0.0;
    }
}


nmMatrix::nmMatrix (int nRows, int mCols, const double * arrayData){
    N_Rows = nRows;
    M_Cols = mCols;
    NM_Elements = N_Rows * M_Cols;
    matrixData = new double [NM_Elements];
    for (int ii = 0; ii < NM_Elements ; ii++){
        *(matrixData + ii) = *(arrayData + ii);
    }
}


nmMatrix::nmMatrix (const nmMatrix & nm_Matrix){
    N_Rows = nm_Matrix.N_Rows;
    M_Cols = nm_Matrix.M_Cols;
    NM_Elements = nm_Matrix.NM_Elements;
    matrixData = new double [NM_Elements];
    for (int ii = 0; ii < NM_Elements ; ii++){
        matrixData[ii] = nm_Matrix.matrixData[ii];
    }
}

nmMatrix::nmMatrix (int nRows, int mCols, char identity){
    N_Rows = nRows; M_Cols = mCols;
    NM_Elements = N_Rows * M_Cols;
    matrixData = new double [NM_Elements];
    if(identity == 'I' || identity == 'i'){
        for(int ii = 0; ii<N_Rows; ii++)
            for (int jj = 0; jj < M_Cols; jj++)
                if (ii == jj)
                    *(matrixData + ii*M_Cols + jj) = 1;
                else *(matrixData + ii*M_Cols + jj) = 0;
    }
}

nmMatrix::~nmMatrix (){
    if(matrixData != nullptr){
        delete [] matrixData;
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%% Elemental Functions %%%%%%%%%%%%%%%%%%%%%%%%%

void nmMatrix::resize (int newRows, int newCols){
    N_Rows = newRows;
    M_Cols = newCols;
    NM_Elements = N_Rows*M_Cols;
    delete[] matrixData;
    matrixData = new double [NM_Elements];
    if (matrixData != nullptr)
        for(int ii = 0 ; ii < NM_Elements ; ii++)
            matrixData[ii] = 0.0;
    else std::cout << "resize nullptr error" << std::endl;
}

double nmMatrix::GetElements (int row, int col){
    int position_e = GetPosition(row, col);
    if ( position_e >= 0){
        return matrixData[position_e];
    }
    else return 0.0/0.0;
}

int nmMatrix::GetPosition (int row, int col){
    if ((row >= 0) && (row < N_Rows) && (col >= 0) && (col < M_Cols)){
        return row*M_Cols + col;
    }
    else return -1;
}

void nmMatrix::replace (int row, int col, const double value){
    int position_e = GetPosition (row, col);
    matrixData [position_e] = value;
}

void nmMatrix::showMatrix(){
    for(int ii = 0 ; ii< this->N_Rows; ii++){
        for(int jj = 0 ; jj< this->M_Cols; jj++){
            std::cout << this->GetElements(ii, jj) << "  ";
        }
        std::cout << "\n";
    }
    std:;cout << "\n";
}

// %%%%%%%%%%%%%%%%%%%%%%%%% Operations on Matrix %%%%%%%%%%%%%%%%%%%%%%%%%

// SUM   %%%%%%%%%%%%%%%%%%%%%%%%%

nmMatrix operator+ (const nmMatrix& lMatrix, const nmMatrix& rMatrix){
    int lRows = lMatrix.N_Rows; int rRows = rMatrix.N_Rows;
    int lCols = lMatrix.M_Cols; int rCols = rMatrix.M_Cols;
    int lElements = lRows*lCols; int rElements = rRows*rCols;
    if((lRows == rRows) && (lCols == rCols)){
        double *partial = new double[lElements];
        for(int ii = 0; ii < lElements; ii++)
            partial[ii] = lMatrix.matrixData[ii] + rMatrix.matrixData[ii];

        nmMatrix result(lRows, lCols, partial);
        delete[] partial;
        return result;
    }
    else {
        std::cout << "Out of range" << std::endl;
        nmMatrix result(1,1);
        return result;
    }
}

nmMatrix operator+= (const nmMatrix& lMatrix, const nmMatrix& rMatrix){
    int lRows = lMatrix.N_Rows; int rRows = rMatrix.N_Rows;
    int lCols = lMatrix.M_Cols; int rCols = rMatrix.M_Cols;
    int lElements = lRows*lCols; int rElements = rRows*rCols;
    if((lRows == rRows) && (lCols == rCols)){
        for(int ii = 0; ii < lElements; ii++)
            lMatrix.matrixData[ii] += rMatrix.matrixData[ii];
        return lMatrix;
    }
    else {
        std::cout << "Out of range" << std::endl;
        nmMatrix result(1,1);
        return result;
    }
}


nmMatrix operator+ (const double & lscaler, const nmMatrix& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows*rCols;
    double *partial_res = new double[rElements];
    for(int ii = 0; ii < rElements; ii++){
        partial_res[ii] = lscaler + rMatrix.matrixData[ii];}

    nmMatrix result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}


nmMatrix operator+ (const nmMatrix& lMatrix, const double & rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    double *partial_res = new double[lElements];
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix.matrixData[ii] + rscaler;}

    nmMatrix result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

nmMatrix operator+= (const nmMatrix& lMatrix, const double & rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    for(int ii = 0; ii < lElements; ii++)
        lMatrix.matrixData[ii] += rscaler;
    return lMatrix;
}

// SUBTRACTION   %%%%%%%%%%%%%%%%%%%%%%%%%


nmMatrix operator- (const nmMatrix& lMatrix, const nmMatrix& rMatrix){
    int lRows = lMatrix.N_Rows; int rRows = rMatrix.N_Rows;
    int lCols = lMatrix.M_Cols; int rCols = rMatrix.M_Cols;
    int lElements = lRows*lCols; int rElements = rRows*rCols;
    if((lRows == rRows) && (lCols == rCols)){
        double *partial = new double[lElements];
        for(int ii = 0; ii < lElements; ii++)
            partial[ii] = lMatrix.matrixData[ii] - rMatrix.matrixData[ii];

        nmMatrix result(lRows, lCols, partial);
        delete[] partial;
        return result;
    }
    else {
        std::cout << "Out of range" << std::endl;
        nmMatrix result(1,1);
        return result;
    }
}

nmMatrix operator-= (const nmMatrix& lMatrix, const nmMatrix& rMatrix){
    int lRows = lMatrix.N_Rows; int rRows = rMatrix.N_Rows;
    int lCols = lMatrix.M_Cols; int rCols = rMatrix.M_Cols;
    int lElements = lRows*lCols; int rElements = rRows*rCols;
    if((lRows == rRows) && (lCols == rCols)){
        for(int ii = 0; ii < lElements; ii++)
            lMatrix.matrixData[ii] -= rMatrix.matrixData[ii];
        return lMatrix;
    }
    else {
        std::cout << "Out of range" << std::endl;
        nmMatrix result(1,1);
        return result;
    }
}

nmMatrix operator- (const double & lscaler, const nmMatrix& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows*rCols;
    double *partial_res = new double[rElements];
    for(int ii = 0; ii < rElements; ii++){
        partial_res[ii] = lscaler - rMatrix.matrixData[ii];}

    nmMatrix result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}


nmMatrix operator- (const nmMatrix& lMatrix, const double & rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    double *partial_res = new double[lElements];
    for(int ii = 0; ii < lElements; ii++){
        partial_res[ii] = lMatrix.matrixData[ii] - rscaler;}

    nmMatrix result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

nmMatrix operator-= (const nmMatrix& lMatrix, const double & rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows*lCols;
    for(int ii = 0; ii < lElements; ii++)
        lMatrix.matrixData[ii] -= rscaler;
    return lMatrix;
}

//MULTIPLICATION (maybe longer forward I could implement multiplication by blocking method)


nmMatrix operator* (const nmMatrix& lMatrix, const nmMatrix& rMatrix){
    int lRows = lMatrix.N_Rows; int lCols = lMatrix.M_Cols;
    int rRows = rMatrix.N_Rows; int rCols = rMatrix.M_Cols;

    if (lCols == rRows){
        double *partial_res = new double[lRows*rCols];
        for(int ii = 0; ii < lRows; ii++){
            for(int jj = 0 ; jj < rCols; jj++){
                double element_res = 0.0;
                for(int kk = 0; kk < lCols; kk++){
                    int linx = ii*lCols + kk;
                    int rinx = kk*rCols + jj;

                    element_res += lMatrix.matrixData[linx] * rMatrix.matrixData[rinx];
                }
                int resx = ii*rCols+jj;

                partial_res[resx] = element_res;
            }
        }

        nmMatrix result(lRows, rCols, partial_res);
        delete[] partial_res;
        return result;
    }
    else {
        std::cout << "Out of range" << std::endl;
        nmMatrix result(1,1);
        return result;
    }
}


nmMatrix operator* (const double & lscaler, const nmMatrix& rMatrix){
    int rRows = rMatrix.N_Rows;
    int rCols = rMatrix.M_Cols;
    int rElements = rRows * rCols;
    double *partial_res = new double[rElements];
    for(int ii = 0 ; ii < rElements; ii++){
        partial_res[ii] = lscaler * rMatrix.matrixData[ii];
    }

    nmMatrix result(rRows, rCols, partial_res);
    delete[] partial_res;
    return result;
}


nmMatrix operator* (const nmMatrix& lMatrix, const double & rscaler){
    int lRows = lMatrix.N_Rows;
    int lCols = lMatrix.M_Cols;
    int lElements = lRows * lCols;
    double *partial_res = new double[lElements];
    for(int ii = 0 ; ii < lElements; ii++){
        partial_res[ii] = lMatrix.matrixData[ii] * rscaler;
    }

    nmMatrix result(lRows, lCols, partial_res);
    delete[] partial_res;
    return result;
}

//EQUALITY   %%%%%%%%%%%%%%%%%%%%%%%%%


bool nmMatrix::operator== (const nmMatrix& nmMatrix){
    if((this->N_Rows != nmMatrix.N_Rows) && (this->M_Cols != nmMatrix.M_Cols))
        return false;

    bool bug = true;
    for(int ii = 0; ii < this->NM_Elements ; ii++){
        if(this->matrixData[ii] != nmMatrix.matrixData[ii]){
            bug = false;
        }
    }

    return bug;
}

// %%%%%%%%%%%%%%%%%%%%%%%%% Special Functions %%%%%%%%%%%%%%%%%%%%%%%%%

void nmMatrix::uniform_random(){
    for(int ii = 0; ii < this->NM_Elements; ii++)
        this->matrixData[ii] = rand64.r();
}

void nmMatrix::normal_random(const double mu, const double sigma){
    for(int ii = 0; ii < this->NM_Elements; ii++)
        this->matrixData[ii] = rand64.gauss(mu, sigma);
}

#endif
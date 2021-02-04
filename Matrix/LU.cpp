#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <vector>
#include <chrono>

using namespace std;

int main()
{
    // int rows = 4;
    // int cols = 4;

    // auto *dense_mat = new Matrix<double>(rows, cols, true);
    // auto *L = new Matrix<double>(rows, cols, true);
    // auto *U = new Matrix<double>(rows, cols, true);

    // vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};

    // for (int i = 0; i < rows * cols; i++)
    // {
    //   dense_mat->values[i] = vs[i];
    // }

    // dense_mat->printMatrix();
    // dense_mat->LUDecomp(*L,*U);
    // L->printMatrix();
    // U->printMatrix();
    
    // delete dense_mat;
    // delete U;
    // delete L;

    //^These are the tests for LUDecomp alone

    int rows = 3;
    int cols = 3;

    auto *dense_mat = new Matrix<double>(rows, cols, true);
    
    vector<double> vs = {4,1,1,4,2,5,0,2,1};
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }

    //auto *LU = new Matrix<double>(rows, cols, true);
    // dense_mat->IPLUDecomp();
    // dense_mat->printMatrix();

    double b[] = {3.0,1.0,3.0};
    double output[] = {0.0,0.0,0.0};

    // // dense_mat->LUSolve(b,output);

    //dense_mat->conjugate_gradient(b,output,100000,0.005);

    dense_mat->LUSolve(b,output,false);
    cout<< "The solution is: ";
    for (int i = 0; i<cols;i++)
    {
      cout<< output[i] << " ";
    }

    // output[0]= 0.0;
    // output[1] = 0.0;
    // output[2] = 0.0;
    
    // dense_mat->LUSolve(b,output, true);

    // cout<< "The solution is: ";
    // for (int i = 0; i<cols;i++)
    // {
    //   cout<< output[i] << " ";
    // }
    // delete dense_mat;
}
#include <iostream>
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <vector>

using namespace std;

int main()
{
    int rows = 4;
    int cols = 4;

    auto *dense_mat = new Matrix<double>(rows, cols, true);
    auto *L = new Matrix<double>(rows, cols, true);
    auto *U = new Matrix<double>(rows, cols, true);

    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};

    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }

    dense_mat->printMatrix();
    dense_mat->LUDecomp(*L,*U);
    L->printMatrix();
    U->printMatrix();
    
    delete dense_mat;
    delete U;
    delete L;
}
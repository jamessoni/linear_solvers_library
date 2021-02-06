#include <iostream>
#include <math.h>
#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

int main_gauss_seidel()
{
    int rows = 0;
    int cols = 0;
    float tol = 1e-6;

    cout << "Input the number of Matrix rows: ";
    cin >> rows;
    //ensuring a square matrix
    cols = rows;

    //creating the SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));

    A->printMatrix();

    //creating space for solution matrix (x) and matrix (b)
    //auto* x_init = new Matrix<double>(rows, 1, true);
    //auto* b = new Matrix<double>(rows, 1, true);

    //// filling x with 0's as initial guess
    //for (int i = 0; i < rows; i++)
    //{
    //    x_init->values[i] = 0;
    //}

    //// filling b with random ints
    //for (int i = 0; i < rows; i++)
    //{
    //    b->values[i] = rand() % 10;
    //}

    double* b = new double[rows];
    for (int i = 0; i < rows; i++) {
        b[i] = rand() % 300 + 5;
    }

    // set up arrays to store solutions and answer check
    double* x_init = new double[rows];
    double* answer_check = new double[rows];
    for (int i = 0; i < rows; i++)
    {
        x_init[i] = 0.0;
    }

    cout << "\nx: ";
    //x_init->printValues();
    //x_init->printMatrix();

    cout << "\nb: ";
    //b->printValues(); 
    //b->printMatrix();

    //using gauss-seidel method
    clock_t start = clock();
    A->gauss_seidel(*A, b, x_init, tol);
    clock_t end = clock();

    cout << "\nSolution x: ";
    //x_init->printValues();
    //x_init->printMatrix();

    cout << "\nTime spent: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " ms" << endl;

    delete[] b;
    delete[] x_init;
    delete[] answer_check;
    delete A;

    return 0;
}
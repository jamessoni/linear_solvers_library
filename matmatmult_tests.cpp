#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <chrono>
#include <memory>
#include <fstream>
#include <string>
#include <sstream> //to split strings
#include <iterator> //to split strings as well
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using namespace std::chrono;

void readMatrixFromFile(string name, Matrix<double> *toread)
{
    vector<string> *splitl;
    fstream myfile;
    string line;
    myfile.open("premade_matrices/"+name);
    while (getline(myfile, line)) {
        istringstream iss(line);
        vector<string> splitline((istream_iterator<string>(iss)), istream_iterator<string>());
        for(int i = 0; i<splitline.size();i++)
        {
            toread->values[i] = stod(splitline[i]);
        }
    }
    myfile.close();
}

vector<int> vs = {10};

void test_matMatMult()
{   
    for (int i=0;i<vs.size();i++)
    {
        cout << "Testing matMatMult():" << "\n";
        int rows = vs[i];
        int cols = vs[i];
        auto *dense_mat1 = new Matrix<double>(rows, cols, true);
        auto *dense_mat2 = new Matrix<double>(rows, cols, true);
        auto *dense_mat3 = new Matrix<double>(rows, cols, true);
        auto *dense_mat4 = new Matrix<double>(rows, cols, true);

        auto *SPDM = new Matrix<double>(3, 3, 10,5);
        SPDM->printMatrix();
        bool spd = SPDM->SPDMatrixcheck();
        std::cout << "\n"<< spd<< "\n";
        SPDM->values[4] = 0;
        SPDM->printMatrix();
        spd = SPDM->SPDMatrixcheck();
        std::cout << "\n" << spd<< "\n";

        readMatrixFromFile("MMM"+ to_string(rows) +"-1.txt",dense_mat1);
        readMatrixFromFile("MMM"+ to_string(rows) +"-2.txt",dense_mat2);
        readMatrixFromFile("MMM"+ to_string(rows) +"-3.txt",dense_mat3);
        
        auto start = high_resolution_clock::now();
        dense_mat1->matMatMult(*dense_mat2,*dense_mat4);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        double RMS = dense_mat1->RMS_norm_diff(dense_mat3->values, dense_mat4->values);
        if (RMS > 1.e-6)
        {
            std::cout << "MatMatMult failed.";
        } else
        {
            std::cout << "MatMatMult successful.";
        }
        std::cout << "Time spent to multiply two " << rows <<"x" << cols<< " matrices: " << duration.count() *0.001 << "ms" << "\n\n";

        delete dense_mat1;
        delete dense_mat2;
        delete dense_mat3;
        delete dense_mat4;
    }
}

int main()
{
    std::cout << "Testing components:" << "\n\n";
    test_matMatMult();
    
    return 0;
}
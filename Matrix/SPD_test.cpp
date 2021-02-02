#include <iostream>
#include <math.h>
#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

//testing both the SPD constructor 
// as well as the SPD (weak) check method
int main()
{
	//testing on small matrices
	int rows = 3;
	int cols = 3;

	//creating the SPD matrix
	auto* A1 = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

	A1->printMatrix();
	A1->SPDMatrixcheck();

	//testing on large matrices
	
	//creating the SPD matrix
	//int rows2 = 1000;
	//int cols2 = 1000;
	//auto* B = new Matrix<double>(rows2, cols2, 3 * pow(rows2, 2), 2 * pow(rows2, 2));

	//B->printMatrix();
	//B->SPDMatrixcheck();

	//delete B;
	delete A1;
	return 0;
}
#include <iostream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <ctime>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"



using namespace std;

int main_22()
{
	vector<double> N;
	vector<double> time_Matrix;
	vector<double> time_element;
	float RMS_a = 0;
	float RMS_b = 0;
	for (int n = 2; n < 1000; n += 10) {

		N.push_back(n);

		cout << "Number of rows and cols " << n << endl;


		int const rows = n;
		int cols = rows;

		// Testing our matrix class
		auto* dense_mat = new Matrix<double>(rows, cols, true);


		// Now let's test printing our matrix with our other function
		//dense_mat->printMatrix();      
		//ddense_mat->printMatrix();



		// create b vector and fill it with random values
		double* foo = new double[rows];
		for (int i = 0; i < rows; i++) {
			foo[i] = rand() % 300 + 5;
		}
		// set up arrays to store solutions for output to screen
		//double output1[rows];
		double* output = new double[rows];
		double* answer_check = new double[n];
		double* answer_check12 = new double[n];



		clock_t start = clock();
		dense_mat->jacobi_solver_matrix(foo, output, 100000, false);
		clock_t end = clock();
		cout << "Time spent on matrix version: " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
		time_Matrix.push_back((double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0);



		dense_mat->matVecMult(output, answer_check12);

		RMS_a = dense_mat->RMS_norm_diff(foo, answer_check12);
		cout << RMS_a << endl;






		clock_t start2 = clock();
		dense_mat->jacobi_solver_element(foo, output, 100000, false);
		clock_t end2 = clock();
		cout << "Time spent element-wise version: " << (double)(end2 - start2) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

		time_element.push_back((double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0);

		// uncomment below check solution, if output to screen is close to zero than last solver is working correct



		dense_mat->matVecMult(output, answer_check12);

		RMS_b = dense_mat->RMS_norm_diff(foo, answer_check12);

		cout << RMS_b << endl;

		delete[] foo;
		delete[] output;
		delete[] answer_check12;
		delete[] answer_check;
		delete[] dense_mat;
	}


	//std::ofstream output_file("./Jacobi_optimized_O3.txt");
	//std::ostream_iterator<std::string> output_iterator(output_file, "\n");
	//std::copy(N.begin(), N.end(), output_iterator);



	system("pause");
	return 0;

}
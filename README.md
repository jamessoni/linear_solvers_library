# Team Void ACSE-5 Group Assignment

### Implementing Linear Solvers

Within this Advanced Programming Assignment the project aimed to implement an algorithm to solve the linear system *A* *x* = *b*, where A is a positive definite matrix, with x and b being vectors.
This linear solver library is written in the programming language, C++, maximising its advantages of Object Oriented Design, memory allocation and optimisation.
The assignment further evaluates, tests and aims to optimise the linear solvers created fulfilling convergence to solutions for both dense and sparse input matrices. Present within this linear solvers library, exists dense and sparse libraries - each with their respective .cpp and .h files. 

#### Getting Started / Pre-requisites
##### Mac/Linux:
IDE or C++ compiler of choice. Recommendations listed below:
* default g++ compiler
* Xcode
##### Windows users:
IDE or C++ compiler of choice. Recommendations listed below:
* Microsoft Visual Studio Community
* Visual Studio Code

##### Installation
Clone or download the [Team Void github repository](https://github.com/acse-2020/group-project-team-void) to your local machine

##### Using
For use of the program will depend on your OS and IDE/compiler of choice:
Visual Studio Community IDE users:
* Open a Visual Studio Project, add all the cloned files to the source files list
* Build project and run.

Notice:
CSRMatrix method matMatMult presents a bug with the use of Visual Studio IDE - stating 'C2760 - syntax error: unexpected token 'identifier', expected ';'.
This bug however is not present using other IDE's. For this reason, we have left the method initially commented.

Terminal:
* Enter directory with cloned repository files within the Matrix/ folder

        $> g++ main_matrix.cpp -fpermissive
        $> .\a
        
The above main_matrix.cpp call will allow the user to enter the User Interface of the Matrix Library. Within such, enables features of the different aspects of the linear solver library to be used. With the users specified input row value, selected solver and on what matrix (dense or sparse); a randomly created SPD matrix and b vector will be used showing the example functionality. Further, you can select the testing framework, which will show a timing comparison of the differing solvers.

#### Linear Solvers:
Within this library there exists various linear solvers; each with their enabled callable method within classes Matrix or CSRMatrix. They can all be accessed and used with the above main_matrix.cpp providing a UI platform to see examples as well as called specifically with the required inputs.

##### Implemented Solvers:
* Jacobi element solver
* Jacobi matrix solver 
* Gauss-seidel solver 
* LU solver
* Conjugate gradient solver
* Cholesky solver

#### Declaring input Matrices:
##### Dense Matrices:
The values of the dense matrices will depend on whether you as the user wish to input a specified matrix.
One example shown below:

    int rows = 4;
    int cols = 4;
    //creating matrix
    auto *dense_mat = new Matrix<double>(rows, cols, true);

    //allocating specific matrix entries
    vector<double> vs = {5, 6, 2, 1, 9, 9, 7, 2, 4, 3, 8, 1, 2, 0, 9, 1};

    //filling Matrix dense_mat with values specified above
    for (int i = 0; i < rows * cols; i++)
    {
      dense_mat->values[i] = vs[i];
    }

    [...]

    delete dense_mat; 
    
The example above declares a dense matrix which has further been filled with the values present in the vector (vs). 
However, in providing ease of use of the system - there is functuality to create a random SPD matrix if necessary.
Shown below:

    int rows = 4;
    int cols = 4;

    //creating SPD matrix
    auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 2), 2 * pow(rows, 2));

    [...]

    delete A;
    
This above with specific argument inputs calls class constructor; in turn producing a randomised Diagonal Dominant Positive Definite Matrix. This constructor call further ensures the matrix produced will converge with specified linear solvers adhering to both diagonal dominance and/or positive definiteness.

##### Sparse Matrices:
Further, sparse matrices like dense share familiar attributes and therefore call from both Matrix and sub-class CSRMatrix.
The sparse matrices are defined through compressed sparse row (CSR) format - pertaining to three One Dimensional arrays that respectively enables fast access of non zero values for methods and calculations.
Example of CSR below:

    int rows = 4;
    int cols = 4;
    int const nnzs = 8;
    
    auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, true);

    int vals[nnzs] = { 10,2,12,2,6,1,1,9 };

    //filling sparse_mat values array
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->values[i] = vals[i];
    }

    //filling sparse_mat col_index array
    int col_ind[nnzs] = { 0,2,1,0,2,3,2,3 };
    for (int i = 0; i < nnzs; i++)
    {
        sparse_mat->col_index[i] = col_ind[i];
    }

    //filling sparse_mat row_position array
    int row_ind[5] = { 0,2,3,6,8 };
    for (int i = 0; i < 5; i++)
    {
        sparse_mat->row_position[i] = row_ind[i];
    }

    [...]
    delete sparse_mat;

##### Solvers
To call the linear solvers (for both dense and sparse) it will require using an object matrix. This will include calling a constructor which randomises its Matrix element values (e.g. SPD constructor) or a constructor where you wish to fill it with values.
Fulfilling the inputs of the solver (refer to the declared solver in header file), you can call a solver as follows:
###### creating inputs SPD Matrix, 1D arrays x and b:

        float tol = 1e-6;
        int rows = 4;
        int cols = 4;
        
        //creating the SPD matrix
        auto* A = new Matrix<double>(rows, cols, 3 * pow(rows, 1.5), 2 * pow(rows, 1.5));

        //creating and filling array b with random values
        double* b = new double[rows];
        for (int i = 0; i < rows; i++) {
            b[i] = rand() % 300 + 5;
        }

        // set up and initialising x (solution) array
        double* x = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            x[i] = 0.0;
        }

###### Example of solver call:
            
        //using gauss-seidel method
        A->gauss_seidel(*A, b, x, tol);
        
#### License
See LICENSE.txt for details.

Testing files, main_tests_components.cpp and main_tests_solvers.cpp are present within the repository for your own leisurely testing.

For further reading/information on specific methods please refer to documentation.md.

###### Team-Void ACSE-5

#include "EigenVectorCalc.h"
#include "Eigen/Eigenvalues"
#include <iostream>
using namespace std;
using namespace Eigen;

void EigenVectorCalc(float **matrix, int matrixsize, float *eigenval, float **eigenvectors){

	MatrixXcf A(matrixsize,matrixsize);
	for (int i = 0; i < matrixsize; ++i)
		for (int j = 0; j < matrixsize; ++j)
			A(i,j) = matrix[i][j];

	ComplexEigenSolver<MatrixXcf> ces;
	ces.compute(A);
	for (int i = 0; i < matrixsize; ++i){
		eigenval[i] = ces.eigenvalues()[i].real();
		for (int j = 0; j < matrixsize; ++j)
		{
			eigenvectors[i][j] = ces.eigenvectors().col(i)[j].real();
		}
	}
}
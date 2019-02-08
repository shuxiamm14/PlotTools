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
	//cout << "The eigenvalues of A are:" << endl << ces.eigenvalues() << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << ces.eigenvectors() << endl << endl;
//
	//cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
	//VectorXcf v = ces.eigenvectors().col(0);
	//cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
	//cout << "... and A * v = " << endl << A * v << endl << endl;
	//cout << "Finally, V * D * V^(-1) = " << endl
	//<< ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << endl;

}
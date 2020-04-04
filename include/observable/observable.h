#ifndef OBSERVABLE
#define OBSERVABLE 1
#include <iostream> 
#include "TH1.h"
class observable {
public: 
    double nominal, error, errordown; 
    observable(double n = 0, double e = 0, double ed = 0);
    void print();

	observable operator = (observable const &obj);
	observable operator + (observable const &obj);
	observable operator - (observable const &obj);
	observable operator += (observable const &obj);
	observable operator -= (observable const &obj);
	observable operator *= (observable const &obj);
	observable operator /= (observable const &obj);
	observable operator * (observable const &obj);
	observable operator / (observable const &obj);
	observable operator + (double aa);
	observable operator - (double aa);
	observable operator * (double aa);
	observable operator / (double aa);

};

observable integral(TH1* histogram, int init = 1, int end = 0);
#endif
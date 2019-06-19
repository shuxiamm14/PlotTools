#include <iostream> 
using namespace std; 
  
class observable { 
public: 
    double nominal, error; 
    observable(double n = 0, double e =0)  {nominal = n;   error = e;} 
      
    void print() { cout << nominal << " +/- " << error << endl; } 

	observable operator = (observable const &obj);
	observable operator + (observable const &obj);
	observable operator - (observable const &obj);
	void operator += (observable const &obj);
	void operator -= (observable const &obj);
	observable operator * (observable const &obj);
	observable operator / (observable const &obj);
	observable operator + (double aa);
	observable operator - (double aa);
	observable operator * (double aa);
	observable operator / (double aa);

};
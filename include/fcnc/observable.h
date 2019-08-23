#include <iostream> 
  
class observable { 
public: 
    double nominal, error; 
    observable(double n = 0, double e =0)  {nominal = n;   error = e;} 
      
    void print() { std::cout << nominal << " +/- " << error << std::endl; } 

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
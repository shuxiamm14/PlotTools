#include "observable.h"
#include "fcnc_include.h"

observable observable::operator + (observable const &obj) { 
	observable res;
	res.nominal = nominal + obj.nominal;
	res.error = rms(error, obj.error);
	res.errordown = res.error;
	return res;
} 

observable observable::operator - (observable const &obj) { 
	observable res;
	res.nominal = nominal - obj.nominal;
	res.error = rms(error, obj.error);
	res.errordown = res.error;
	return res;
} 
observable observable::operator * (observable const &obj) { 
	observable res;
	res.nominal = nominal * obj.nominal;
	res.error = rms(error * obj.nominal, obj.error * nominal);
	res.error = rms(res.error, error*obj.error);
	res.errordown = res.error;
	return res;
} 
observable observable::operator / (observable const &obj) { 
	observable res;
	res.nominal = nominal / obj.nominal;
	res.error = rms(error / obj.nominal, obj.error * nominal / obj.nominal / obj.nominal);
	res.errordown = res.error;
	return res;
} 

observable observable::operator = (observable const &obj) { 
	nominal = obj.nominal;
	error = obj.error;
	errordown = obj.errordown;
	return *this;
} 

observable observable::operator += (observable const &obj) { 
	*this = *this + obj;
	return *this;
} 

observable observable::operator -= (observable const &obj) { 
	*this = *this - obj;
	return *this;
} 

observable observable::operator *= (observable const &obj) { 
	*this = *this * obj;
	return *this;
} 

observable observable::operator /= (observable const &obj) { 
	*this = *this / obj;
	return *this;
} 

observable observable::operator + (double aa) { 
	observable res;
	res.nominal = nominal + aa;
	res.error = error;
	res.errordown = errordown;
	return res;
} 
observable observable::operator - (double aa) { 
	observable res;
	res.nominal = nominal - aa;
	res.error = error;
	res.errordown = errordown;
	return res;
} 
observable observable::operator * (double aa) { 
	observable res;
	res.nominal = nominal * aa;
	res.error = error * aa;
	res.errordown = errordown * aa;
	return res;
} 
observable observable::operator / (double aa) { 
	observable res;
	res.nominal = nominal / aa;
	res.error = error / aa;
	res.errordown = errordown / aa;
	return res;
} 

observable integral(TH1* histogram, int init, int end)
{
	if(end == 0) end = histogram->GetNbinsX();
	double err;
	double itg(histogram->IntegralAndError(init, end, err));
	observable ret(itg,err);
	return ret;
}

observable::observable(double n, double e, double ed) : nominal(n), error(e) {
	if(ed == 0) errordown = e;
	else errordown = ed;
}
void observable::print(){
	printf("%f+%f-%f",nominal,error,errordown);
}

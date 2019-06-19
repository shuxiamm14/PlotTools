#include "observable.h"
#include "fcnc_include.h"

observable observable::operator + (observable const &obj) { 
     observable res; 
     res.nominal = nominal + obj.nominal; 
     res.error = rms(error, obj.error); 
     return res; 
} 

observable observable::operator - (observable const &obj) { 
     observable res; 
     res.nominal = nominal - obj.nominal; 
     res.error = rms(error, obj.error); 
     return res; 
} 
observable observable::operator * (observable const &obj) { 
     observable res; 
     res.nominal = nominal * obj.nominal; 
     res.error = rms(error * obj.nominal, obj.error * nominal); 
     return res; 
} 
observable observable::operator / (observable const &obj) { 
     observable res; 
     res.nominal = nominal / obj.nominal; 
     res.error = rms(error / obj.nominal, obj.error * nominal / obj.nominal / obj.nominal); 
     return res; 
} 

observable observable::operator = (observable const &obj) { 
     nominal = obj.nominal; 
     error = obj.error; 
     return *this; 
} 

void observable::operator += (observable const &obj) { 
     nominal = nominal + obj.nominal; 
     error = rms(error, obj.error); 
} 

void observable::operator -= (observable const &obj) { 
     nominal = nominal - obj.nominal; 
     error = rms(error, obj.error); 
} 

observable observable::operator + (double aa) { 
     observable res; 
     res.nominal = nominal + aa; 
     res.error = error; 
     return res; 
} 
observable observable::operator - (double aa) { 
     observable res; 
     res.nominal = nominal - aa; 
     res.error = error; 
     return res; 
} 
observable observable::operator * (double aa) { 
     observable res; 
     res.nominal = nominal * aa; 
     res.error = error * aa; 
     return res; 
} 
observable observable::operator / (double aa) { 
     observable res; 
     res.nominal = nominal / aa; 
     res.error = error / aa; 
     return res; 
} 
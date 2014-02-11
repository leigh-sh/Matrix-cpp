#include <cmath>
#include "ExNum.h"

ExNum::ExNum(double source){
	_value = source;
	_value *= MUL_FACTOR;
	if(_value > 0){
		_value = floor(_value);
	}
	else{
		_value = ceil(_value);
	}
	_value /= MUL_FACTOR;
}

ExNum&
ExNum::operator=(const ExNum &other){
	if(&other == this){
		return *this;
	}
	_value = other._value;
	return *this;
}

ExNum
const ExNum::operator+(const ExNum &other) const{
	ExNum res(_value + other._value);
	return res;
}

ExNum
const ExNum::operator-(const ExNum &other) const{
	ExNum res(_value - other._value);
	return res;
}

ExNum
const ExNum::operator*(const ExNum &other) const{
	ExNum res(_value * other._value);
	return res;
}

ExNum
const ExNum::operator/(const ExNum &other) const{
	ExNum res(_value / other._value);
	return res;
}

ExNum
const ExNum::operator-() const{
	if(_value != 0){
		ExNum temp(-_value);
		return temp;
	}
	return (*this);
}

ExNum&
ExNum::operator+=(const ExNum &other){
	*this = *this + other;
	return *this;
}

ExNum&
ExNum::operator-=(const ExNum &other){
	*this = *this - other;
	return *this;
}

ExNum&
ExNum::operator*=(const ExNum &other){
	*this = *this * other;
	return *this;
}

ExNum&
ExNum::operator/=(const ExNum &other){
	*this = *this / other;
	return *this;
}
bool 
ExNum::operator==(const ExNum &other) const{
	return (_value == other._value);
}

bool 
ExNum::operator!=(const ExNum &other) const{
	return (_value != other._value);
}

double 
ExNum::getValue() const{
	return _value;
}

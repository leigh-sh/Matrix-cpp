#ifndef EXNUM_H
#define EXNUM_H

#define MUL_FACTOR 100

/**
* This class is a double of only 2 digits after the floating point, with
* all relevant operators.
*/

class ExNum{

public:
	ExNum(double source);
	ExNum& operator=(const ExNum& other);
	const ExNum operator+(const ExNum &other) const;
	const ExNum operator-(const ExNum &other) const;
	const ExNum operator*(const ExNum &other) const;
	const ExNum operator/(const ExNum &other) const;
	const ExNum operator-() const;
	ExNum& operator+=(const ExNum &other);
	ExNum& operator-=(const ExNum &other);
	ExNum& operator*=(const ExNum &other);
	ExNum& operator/=(const ExNum &other);
	bool operator==(const ExNum &other) const;
	bool operator!=(const ExNum &other) const;
	double getValue() const;

private:
	double _value;
};

#endif

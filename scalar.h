#ifndef SCALAR_H
#define SCALAR_H

#include "ExNum.h"

/**
* This class describes a Scalar - a coordinate in a sparse matrix, with 
* members - row, col and the value itself (of type ExNum).
*/

class Scalar{

public:
	Scalar(double value, int row, int col);
	Scalar(ExNum num, int row, int col);
	ExNum getNum() const;
	int getRow() const;
	int getCol() const;
	void setRow(int row);
	void setCol(int col);
	const Scalar operator+(const Scalar& other) const;
	const Scalar operator-() const;
	bool operator==(const Scalar &other) const;
	bool operator!=(const Scalar &other) const;
	void transpose();

private:
	ExNum _num;
	int _row;
	int _col;
};

#endif

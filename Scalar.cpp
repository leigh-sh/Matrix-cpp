#include <cassert>
#include "Scalar.h"

Scalar::Scalar(double value, int row, int col):
	_num(value), _row(row), _col(col)
{

}

Scalar::Scalar(ExNum num, int row, int col):
	_num(num), _row(row), _col(col)
{

}

ExNum
Scalar::getNum() const{
	return _num;
}

int
Scalar::getRow() const{
		return _row;
}

int
Scalar::getCol() const{
	return _col;
}

void 
Scalar::setRow(int row){
	_row = row;
}

void 
Scalar::setCol(int col){
	_col = col;
}

Scalar
const Scalar::operator+(const Scalar& other) const{
	assert(_row == other._row && _col == other._col);
	Scalar temp(_num.getValue() + other._num.getValue(), _row, _col);
	return temp;
}

Scalar
const Scalar::operator-() const{
	if(_num == ExNum(0)){
		return (*this);
	}
	Scalar temp(-_num.getValue(), _row, _col);
	return temp;
}

bool
Scalar::operator==(const Scalar &other) const{
	return (_row == other._row &&
			_col == other._col &&
			_num == other._num);
}

bool
Scalar::operator!=(const Scalar &other) const{
	return !((*this) == other);
}

void
Scalar::transpose(){
	int temp = _row;
	_row = _col;
	_col = temp;
}

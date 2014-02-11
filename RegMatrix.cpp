#include <iostream>
#include <cassert>
#include <cmath>
#include "ExNum.h"
#include "RegMatrix.h"
#include "SparseMatrix.h"
#include "Scalar.h"

//Constructor. if fill == NULL, creates a RegMatrix of given dimensions
//with only zeroes
RegMatrix::RegMatrix(int numRow, int numCol, double fill[]):
	_numRow(numRow), _numCol(numCol)
{
	assert(_numRow > 0 && _numCol > 0);
	if(NULL != fill){
		for(int i=0; i<numRow; ++i){
			std::vector<ExNum> v;
			for(int j=0; j<numCol; ++j){
				ExNum num(fill[i*numCol + j]);
				v.push_back(num);
			}
			_vec.push_back(v);
		}
	}
	else{
		(*this) = RegMatrix(_numRow, _numCol);
	}
}

//Constructor as described. Assumption: &sparse != NULL
RegMatrix::RegMatrix(const SparseMatrix &sparse)
{
	assert(NULL != &sparse);
	*this = RegMatrix(sparse.getNumRow(), sparse.getNumCol());
	//Getting a RegMatrix of the right dimensions, with only zeroes.
	assert(_numRow > 0 && _numCol > 0);
	std::vector<Scalar> v = sparse.getVector();
	for(unsigned int i=0; i<v.size(); ++i){//filling the RegMatrix
		_vec[v[i].getRow()][v[i].getCol()] = v[i].getNum();
	}
}

//A constructor creating a RegMatrix of the given dimensions with only zeroes
RegMatrix::RegMatrix(int numRow, int numCol):
	_numRow(numRow), _numCol(numCol)
{
	assert(_numRow > 0 && _numCol > 0);
	for(int i=0; i<numRow; ++i){
		std::vector<ExNum> v;
		for(int j=0; j<numCol; ++j){
			ExNum num(0);
			v.push_back(num);
		}
		_vec.push_back(v);
	}
}

//Returning the minor of this matrix, without the given row and col
RegMatrix
RegMatrix::minor(int delRow, int delCol) const{
	RegMatrix minor(_numRow - 1, _numCol - 1);
	for(int i=0; i<minor._numRow; ++i){
		for(int j=0; j<minor._numCol; ++j){
			int i2, j2;
			if(i < delRow){
				i2 = i;
			}
			else{
				i2 = i+1;
			}
			if(j < delCol){
				j2 = j;
			}
			else{
				j2 = j + 1;
			}
			minor._vec[i][j] = _vec[i2][j2];
		}
	}
	return minor;
}

int 
RegMatrix::getNumRow() const{
	return _numRow;
}

int 
RegMatrix::getNumCol() const{
	return _numCol;
}

//Returns the num in the given coordinates
ExNum
RegMatrix::numAt(int row, int col) const{
	return _vec[row][col];
}

//Sets the number in the given coordinates
void 
RegMatrix::setNumAt(int row, int col, ExNum num){
	_vec[row][col] = num;
}

RegMatrix
const RegMatrix::operator+(const RegMatrix& other) const{
	assert(_numRow == other._numRow && _numCol == other._numCol);
	RegMatrix temp(_numRow, _numCol);
	for(int i=0; i<_numRow; ++i){
		for(int j=0; j<_numCol; ++j){
			temp._vec[i][j] = _vec[i][j] + other._vec[i][j];
		}
	}
	return temp;
}

RegMatrix
const RegMatrix::operator*(const RegMatrix& other) const{
	assert(_numCol == other._numRow);
	RegMatrix temp(_numRow, other._numCol);
	for(int i=0; i<temp._numRow; ++i){
		for(int j=0; j<temp._numCol; ++j){
			for(int h=0; h<_numCol; ++h){
				ExNum n = _vec[i][h] * other._vec[h][j];
				temp._vec[i][j] += n;
			}
		}
	}
	return temp;
}

RegMatrix
const RegMatrix::operator-() const{
	RegMatrix temp(_numRow, _numCol);
	for(int i=0; i<temp._numRow; ++i){
		for(int j=0; j<temp._numCol; ++j){
			temp._vec[i][j] = -_vec[i][j];
		}
	}
	return temp;
}

RegMatrix& 
RegMatrix::operator+=(const RegMatrix& other){
	*this = *this + other;
	return *this;
}

RegMatrix& 
RegMatrix::operator+=(const SparseMatrix& sparse){
	*this = *this + sparse;
	return *this;
}

RegMatrix& 
RegMatrix::operator*=(const RegMatrix& other){
	*this = *this * other;
	return *this;
}

RegMatrix& 
RegMatrix::operator*=(const SparseMatrix& sparse){
	*this = *this * sparse;
	return *this;
}

double 
RegMatrix::sum() const{
	ExNum n(0);
	for(int i=0; i<_numRow; ++i){
		for(int j=0; j<_numCol; ++j){
			n += _vec[i][j];
		}
	}
	return n.getValue();
}

bool 
RegMatrix::operator==(const RegMatrix &other) const{
	assert(_numRow == other._numRow && _numCol == other._numCol);
	for(int i=0; i<_numRow; ++i){
		for(int j=0; j<_numCol; ++j){
			if(_vec[i][j] != other._vec[i][j]){
				return false;
			}
		}
	}
	return true;
}

bool 
RegMatrix::operator!=(const RegMatrix &other) const{
	return !(*this == other);
}

RegMatrix 
RegMatrix::transpose() const{
	RegMatrix temp(_numCol, _numRow);
	for(int i=0; i<temp._numRow; ++i){
		for(int j=0; j< temp._numCol; ++j){
			temp._vec[i][j] = _vec[j][i];
		}
	}
	return temp;
}

//Calculates the RegMatrix det according to the first row, recursively
bool
RegMatrix::det(double& result){
	if(_numRow != _numCol){
		return false;
	}
	if(1 == _numRow || 2 == _numRow){
		result = baseDet();
		return true;
	}
	else{
		ExNum n(0);
		for(int j=0; j<_numCol; ++j){
			RegMatrix minor = this->minor(0, j);
			double detMinor = 0;
			minor.det(detMinor);
			n += ExNum(pow(-1.0, j + 2)) * ExNum(_vec[0][j]) * ExNum(detMinor);
		}
		result = n.getValue();
		return true;
	}
}

//Returns the det of this RegMatrix iff it is of size 1*1 or 2*2,
//Therefore this method is of const complexity
double 
RegMatrix::baseDet(){
	assert(_numRow == _numCol);
	assert(_numRow == 1 || _numRow == 2); 
	if(1 == _numRow){
		return _vec[0][0].getValue();
	}
	if(2 == _numRow){
		return ((_vec[0][0] * _vec[1][1]) - (_vec[0][1] * _vec[1][0])).getValue();
	}
        return 0;
}

//Prints the RegMatrix
void
RegMatrix::print() const{
    //	std::cout << "A RegMatrix of dimensions "<< _numRow << "X" << _numCol << "\n";
	for(unsigned int i=0; i<_vec.size(); ++i){
		for(unsigned int j=0; j<_vec[i].size(); ++j){
			std::cout << _vec[i][j].getValue() << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

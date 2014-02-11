#include <iostream>
#include <cassert>
#include <cmath>
#include "SparseMatrix.h"

//Constructor. if fill == NULL, creates a SparseMatrix of given dimensions
//with only zeroes
SparseMatrix::SparseMatrix(int numRow, int numCol, double fill[]):
	_numRow(numRow), _numCol(numCol)
{
	assert(_numRow > 0 && _numCol > 0);
	if(NULL != fill){
		for(int i=0; i<numRow; ++i){
			for(int j=0; j<numCol; ++j){
				if(ExNum(0) != ExNum(fill[i*numCol + j])){
					Scalar sc(fill[i*numCol + j], i, j);
					_vec.push_back(sc);
				}
			}
		}
	}
}

//Constructor as described. Assumption: regMat != NULL
SparseMatrix::SparseMatrix(const RegMatrix &regMat):
	_numRow(regMat.getNumRow()), _numCol(regMat.getNumCol())
{
	assert(_numRow > 0 && _numCol > 0);
	for(int i=0; i<_numRow; ++i){
		for(int j=0; j<_numCol; ++j){
			ExNum num = regMat.numAt(i, j);
			if(ExNum(0) != num){
				Scalar sc(num, i, j);
				_vec.push_back(sc);
			}
		}
	}
}

//Returning the minor of this matrix, without the given row and col
SparseMatrix 
SparseMatrix::minor(int delRow, int delCol) const{
	assert(delRow < _numRow && delCol < _numCol);
	SparseMatrix minor = (*this);
	minor._numRow = _numRow - 1;
	minor._numCol = _numCol - 1;
	for(int i=0; i<minor.size(); ++i){
		if(minor._vec[i].getRow() == delRow || minor._vec[i].getCol() == delCol){
			minor._vec.erase(minor._vec.begin() + i);
			--i;
		}
	}
	for(int i=0; i<minor.size(); ++i){
		if(minor._vec[i].getRow() > delRow){
			minor._vec[i].setRow(minor._vec[i].getRow() - 1);
		}
		if(minor._vec[i].getCol() > delCol){
			minor._vec[i].setCol(minor._vec[i].getCol() - 1);
		}
	}
	return minor;
}

int
SparseMatrix::getNumRow() const{
	return _numRow;
}

int
SparseMatrix::getNumCol() const{
	return _numCol;
}

int 
SparseMatrix::size() const{
	return _vec.size();
}

//Returns a vector with scalars of this SparseMatrix
std::vector<Scalar>
SparseMatrix::getVector() const{
	return _vec;
}

double 
SparseMatrix::sum() const{
	ExNum n(0);
	for(int i=0; i<size(); ++i){
		n += _vec[i].getNum();
	}
	return n.getValue();
}

SparseMatrix
const SparseMatrix::operator+(const SparseMatrix& other) const{
	assert(_numRow == other._numRow && _numCol == _numCol);
	SparseMatrix temp(_numRow, _numCol, NULL);
	int i=0, j=0;
	while(i < size() && j<other.size()){
		if(_vec[i].getRow() < other._vec[j].getRow() ||
			(_vec[i].getRow() == other._vec[j].getRow() && _vec[i].getCol() < other._vec[j].getCol())){
				temp._vec.push_back(_vec[i]);
				++i;
				//The coordinate in this SparseMatrix precedes the one in other
		}
		else if(_vec[i].getRow() == other._vec[j].getRow() && _vec[i].getCol() == other._vec[j].getCol()){
			Scalar sc((_vec[i].getNum() + other._vec[j].getNum()).getValue(), _vec[i].getRow(), _vec[i].getCol());
			if(sc.getNum() != ExNum(0)){
				temp._vec.push_back(sc);
				//The specific coordinate in both matrices are not 0, and the sum of both is not zero
			}
			++i;
			++j;
		}
		else{
			temp._vec.push_back(other._vec[j]);
			++j;
			//The coordinate in other SparseMatrix precedes the one in this
		}
	}
	//Naturally maximum 1 of the next two loops will take place!
	while(i < size()){
		Scalar sc(_vec[i].getNum().getValue(), _vec[i].getRow(), _vec[i].getCol());
		temp._vec.push_back(sc);
		++i;
	}
	while(j < other.size()){
		Scalar sc(other._vec[j].getNum().getValue(), other._vec[j].getRow(), other._vec[j].getCol());
		temp._vec.push_back(sc);
		++j;
	}
	return temp;
}

//Returns true iff the value of [row][col] in this sparseMatrix is not 0,
//and changes index to be the index of it in vec. otherWise index will be 
//the following coordinate
//Assumption: the Scalar[row][col] will come only after(!) index in the 
//order of vec (saves complexity etc).
bool 
SparseMatrix::find(int row, int col, int &index){
	while(index < size() && 
		((_vec[index].getRow() < row) || (_vec[index].getRow() == row && _vec[index].getCol() < col))){
			++index;
	}
	if(index == size()){
		return false;
	}
	if(_vec[index].getRow() == row && _vec[index].getCol() == col){
		return true;
	}
	return false;
}

SparseMatrix
const SparseMatrix::operator*(const SparseMatrix& other) const{
	assert(_numCol == other._numRow);
	SparseMatrix temp(_numRow, other._numCol, NULL);
	SparseMatrix trans = other.transpose();
	int index1 = 0, index2 = 0;
	//The easy way to multiply: transposing the second SparseMatrix
	//and multiplying row by row instead of row by col
	for(int i=0; i<temp._numRow; ++i){
		while(index1 < size() && _vec[index1].getRow() < i){
			++index1;//pushing until we are in the correct row
		}
		for(int j=0; j<temp._numCol; ++j){
			ExNum n(0);
			int startRow = index1;
			while(index1 < size() &&_vec[index1].getRow() == i){
				if(trans.find(j, _vec[index1].getCol(), index2)){
					//finding the coordinate in trans in which we should multiply
					n += _vec[index1].getNum() * trans._vec[index2].getNum();
					++index2;
				}
				++index1;
			}
			if(n!= ExNum(0)){
				temp._vec.push_back(Scalar(n, i, j));
			}
			index1 = startRow;
		}
		index2 = 0;
	}
	return temp;
}

SparseMatrix
const SparseMatrix::operator-() const{
	SparseMatrix temp = (*this);
	for(int i=0; i<size(); ++i){
		temp._vec[i] = -_vec[i];
	}
	return temp;
}

SparseMatrix&
SparseMatrix::operator+=(const SparseMatrix& other){
	(*this) = (*this) + other;
	return (*this);
}

SparseMatrix& 
SparseMatrix::operator+=(const RegMatrix& reg){
	*this = SparseMatrix(RegMatrix(*this) + reg);
	return *this;
}

SparseMatrix& 
SparseMatrix::operator*=(const SparseMatrix& other){
	(*this) = (*this) * other;
	return (*this);
}

SparseMatrix& 
SparseMatrix::operator*=(const RegMatrix& reg){
	*this = SparseMatrix(RegMatrix(*this) * reg);
	return *this;
}

//Assumption:both matrices' vecs are sorted by rows and cols
bool
SparseMatrix::operator==(const SparseMatrix &other) const{
	assert(_numRow == other._numRow && _numCol == other._numCol);
	if(size() != other.size()){
		return false;
	}
	for(int i=0; i<size(); ++i){
		if(_vec[i] != other._vec[i]){
			return false;
		}
	}
	return true;
}

bool
SparseMatrix::operator!=(const SparseMatrix &other) const{
	return !((*this) == other);
}

RegMatrix operator+(const RegMatrix &reg, const SparseMatrix &sparse){
	assert(reg.getNumRow() == sparse._numRow && reg.getNumCol() == sparse._numCol);
	return reg + RegMatrix(sparse);
}

RegMatrix operator+(const SparseMatrix &sparse, const RegMatrix &reg){
	return (reg + sparse);
}

RegMatrix operator*(const RegMatrix &reg, const SparseMatrix &sparse){
	assert(reg.getNumCol() == sparse._numRow);
	return reg * RegMatrix(sparse);
}

RegMatrix operator*(const SparseMatrix &sparse, const RegMatrix &reg){
	assert(sparse._numCol == reg.getNumRow());
	return RegMatrix(sparse) * reg;
}

SparseMatrix 
SparseMatrix::transpose() const{
	SparseMatrix temp = (*this);
	int t = temp._numRow;//swap rows and cols
	temp._numRow = temp._numCol;
	temp._numCol = t;
	for(int i=0; i<temp.size(); ++i){
		temp._vec[i].transpose();
	}
	temp.sortVec(0);
	return temp;
}

//Sorting the scalars to correct order after transposing
void
SparseMatrix::sortVec(int start){
	if(start == size()){
		return;
	}
	int min = start;
	for(int i=start+1; i<size(); ++i){
		if(_vec[i].getRow() < _vec[min].getRow() || 
			(_vec[i].getRow() == _vec[min].getRow() && _vec[i].getCol() < _vec[min].getCol())){
			min = i;
		}
	}
	Scalar temp = _vec[start];
	_vec[start] = _vec[min];
	_vec[min] = temp;
	sortVec(start + 1);
}

//Calculates the SparseMatrix det according to the first row, recursively
bool
SparseMatrix::det(double& result){
	if(_numRow != _numCol){
		return false;
	}
	if(1 == _numRow || 2 == _numRow){
		result = baseDet();
		return true;
	}
	else{
		ExNum n(0);
		int i=0;
		while(i < size() && _vec[i].getRow() == 0){
			SparseMatrix minor = this->minor(0, _vec[i].getCol());
			double detMinor = 0;
			minor.det(detMinor);
			n += ExNum(pow(-1.0, _vec[i].getCol() + 2)) * _vec[i].getNum() * ExNum(detMinor);
			++i;
		}
		result = n.getValue();
		return true;
	}
}

//Returns the det of this SparseMatrix iff it is of size 1*1 or 2*2,
//Therefore this method is of const complexity
double 
SparseMatrix::baseDet() const{
	assert(_numRow == _numCol);
	assert(1 == _numRow || 2 == _numRow);
	if(1 == _numRow){
		if(0 == size()){
			return 0;
		}
		else{
			return _vec[0].getNum().getValue();
		}
	}
	std::vector<ExNum> temp(4,ExNum(0));
	for(int i=0; i<size(); ++i){
		temp[_vec[i].getRow() * 2 + _vec[i].getCol()] = _vec[i].getNum();
	}
	return ((temp[0] * temp[3]) - (temp[1] * temp[2])).getValue();
}

//Prints the SparseMatrix
void
SparseMatrix::print() const{
	int h = 0;
        //	std::cout << "A SparseMatrix of dimensions "<< _numRow << "X" << _numCol << "\n" <<
	//				"Containing " << size() <<" terms different than 0\n";
	for(int i=0; i<_numRow; ++i){
		for(int j=0; j<_numCol; ++j){
			if(h<size()){
				if(_vec[h].getRow() == i && _vec[h].getCol() == j){
					assert(_vec[h].getNum() != ExNum(0));
					std::cout << _vec[h].getNum().getValue() << " ";
					++h;
				}
				else{
					std::cout << 0 << " ";
				}
			}
			else{
				std::cout << 0 << " ";
			}
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

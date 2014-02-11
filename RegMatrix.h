#ifndef REGMATRIX_H
#define REGMATRIX_H

#include <vector>
#include "ExNum.h"

/**
* A class of RegMatrix, as described in the exercise description.
* Note: rows and cols start from 0, not from 1.
*/

typedef std::vector<ExNum> coolVector;

class SparseMatrix;

class RegMatrix{

public:
	RegMatrix(int numRow, int numCol, double fill[]);
	RegMatrix(const SparseMatrix &sparse);
	RegMatrix(int numRow, int numCol);
	RegMatrix minor(int delRow, int delCol) const;
	int getNumRow() const;
	int getNumCol() const;
	ExNum numAt(int row, int col) const;
	void setNumAt(int row, int col, ExNum num);
	const RegMatrix operator+(const RegMatrix& other) const;
	const RegMatrix operator*(const RegMatrix& other) const;
	const RegMatrix operator-() const;
	RegMatrix& operator+=(const RegMatrix& other);
	RegMatrix& operator+=(const SparseMatrix& sparse);
	RegMatrix& operator*=(const SparseMatrix& sparse);
	RegMatrix& operator*=(const RegMatrix& other);
	bool operator==(const RegMatrix &other) const;
	bool operator!=(const RegMatrix &other) const;
	double sum() const;
	RegMatrix transpose() const;
	bool det(double& result);


	void print() const;

private:
	double baseDet();
	std::vector<coolVector> _vec;
	//A vector holding the matrix coordinates data
	int _numRow;
	int _numCol;
};

#endif

#ifndef SPARSE_H
#define SPARSE_H

#include <vector>
#include "Scalar.h"
#include "RegMatrix.h"

/**
* A class of SparseMatrix, as described in the exercise description.
* Note: rows and cols start from 0, not from 1.
*/

class SparseMatrix{
	
public:
	SparseMatrix(int numRow, int numCol, double fill[]);
	SparseMatrix(const RegMatrix &regMat);
	SparseMatrix minor(int delRow, int delCol) const;
	int getNumRow() const;
	int getNumCol() const;
	int size() const;
	std::vector<Scalar> getVector() const;
	const SparseMatrix operator+(const SparseMatrix& other) const;
	const SparseMatrix operator*(const SparseMatrix& other) const;
	const SparseMatrix operator-() const;
	SparseMatrix& operator+=(const SparseMatrix& other);
	SparseMatrix& operator+=(const RegMatrix& reg);
	SparseMatrix& operator*=(const RegMatrix& reg);
	SparseMatrix& operator*=(const SparseMatrix& other);
	bool operator==(const SparseMatrix &other) const;
	bool operator!=(const SparseMatrix &other) const;
	friend RegMatrix operator+(const RegMatrix &reg, const SparseMatrix &sparse);
	friend RegMatrix operator+(const SparseMatrix &sparse, const RegMatrix &reg);
	friend RegMatrix operator*(const RegMatrix &reg, const SparseMatrix &sparse);
	friend RegMatrix operator*(const SparseMatrix &sparse, const RegMatrix &reg);
	double sum() const;
	SparseMatrix transpose() const;
	bool det(double& result);
	void print() const;

private:
	void sortVec(int start);
	double baseDet() const;
	bool find(int row, int col, int &index);
	std::vector<Scalar> _vec;
	int _numRow;
	int _numCol;
};

#endif

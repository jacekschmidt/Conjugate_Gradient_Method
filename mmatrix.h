#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>

// Class that represents a mathematical matrix
class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};

inline MVector operator*(const MMatrix& A, const MVector& x)
{
	if (x.size()!=A.Cols())
	{
		std::cout<<"incompatible operator"<<std::endl;
		return {0};
	}
	MVector v(A.Rows());
	for (int i=0; i<A.Rows(); i++)
	{
		v[i]=0;
		for(int j=0; j<A.Cols();j++)
		{
			v[i]+=A(i,j)*x[j];
		}
	}
	return v;
}

std::ostream& operator<<(std::ostream& os, const MMatrix& A)
{
	os.width(1);
	os << "(";
	os.width(9);
	os.precision(9);
	os<<A(0,0);
	os<<",";

	for(int i=0;i<A.Rows();i++)
	{
		for(int j=0;j<A.Cols();j++)
		{
			if (i==0 && j==0) continue;
			os.width(10);
			os.precision(10);
			os<<A(i,j);
			if(j<A.Cols()-1) os<<",";
		}
		if(i<A.Rows()-1) os<<std::endl;
	}
	os << ")";
	return os;
}

#endif

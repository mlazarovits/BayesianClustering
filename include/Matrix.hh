#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include "PointCollection.hh"
using std::vector;

class Matrix{
	public:
		//set dimensionality
		Matrix();
		//x_dim = dim (rows), y_dim = # points (cols)
		Matrix(int row, int col);
		Matrix(vector<double> in);
		//copy constructor
		Matrix(const Matrix& mat);
		virtual ~Matrix();
		void InitRandom(double min = 0, double max = 1., unsigned long long seed = 123);
		void InitRandomSym(double min = 0, double max = 1., unsigned long long seed = 123);
		void InitRandomSymPos(double min = 0, double max = 1., unsigned long long seed = 123);
		void InitRandomSymPosDef(double min = 0, double max = 1., unsigned long long seed = 333);
		void InitEmpty();
		//initializes identity matrix
		void InitIdentity();
		//fills this matrix with data pts from n-dim gaus
		void SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample);
		void SetDims(int row, int col);
		void SetEntry(double val, int i, int j);
		double at(int i, int j) const;
		void GetCofactor(Matrix& temp, int p, int q, int n) const;
		vector<int> GetDims() const;
		//multiplies mat1*mat2 and stores result in this
		void mult(const Matrix& mat1, const Matrix& mat2);
		//multiplies mat by factor and stores in this
		void mult(const Matrix& mat, double factor);
		//add two matrices and store result in this
		void add(const Matrix& mat1, const Matrix& mat2);
		//add mat to this
		void add(const Matrix& mat);
		//double det() const;
		double det(int n = 0) const;
		bool square() const { return _square; }
		bool symmetric() const;
		void transpose(const Matrix& mat);
		void adjoint(const Matrix& mat);
		void invert(const Matrix& mat);
		double trace();
		Matrix cholesky();
		void eigenCalc(vector<double>& vals, vector<Matrix>& vecs);
		PointCollection MatToPoints();
		void PointsToMat(PointCollection& pc);
		//clear
		void clear();
		void Print() const;
		void reset(){ clear(); InitEmpty(); }
		//update entries
		//add entries
	private:
		int m_row;
		int m_col;
		vector<vector<double>> m_entries;
		bool _square;
};


#endif

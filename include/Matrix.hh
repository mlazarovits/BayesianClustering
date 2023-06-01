#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include "PointCollection.hh"
using std::vector;

class Matrix{
	public:
		//set dimensionality
		Matrix();
		//x_dim = # pts (rows), y_dim = dimensionality (cols)
		Matrix(int row, int col);
		Matrix(vector<double> in);
		//constructor from inputs
		virtual ~Matrix();

		void InitRandom(double min = 0, double max = 0);	
		void InitRandomSym(double min = 0, double max = 1.);
		void InitRandomSymPosDef(double min = 0, double max = 1.);
		void InitEmpty();
		//fills this matrix with data pts from n-dim gaus
		void SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample);
		void SetDims(int row, int col);
		void SetEntry(double val, int i, int j);
		double at(int i, int j) const;
		void GetCofactor(Matrix& temp, int p, int q, int n) const;
		vector<int> GetDims() const;
		void mult(const Matrix& mat1, const Matrix& mat2);
		//add two matrices and store result in this
		void add(const Matrix& mat1, const Matrix& mat2);
		//add mat to this
		void add(const Matrix& mat);
		//double det() const;
		double det(int n = 0) const;
		bool square() const { return _square; }
		bool symmetric() const;
		void transpose(const Matrix& mat);
		void adjoint(Matrix& adj);
		void invert(Matrix& inv);
		Matrix cholesky();
		PointCollection MatToPoints();
		//clear
		void clear();
		void print() const;
		//update entries
		//add entries
	private:
		int m_row;
		int m_col;
		vector<vector<double>> m_entries;
		bool _square;
};


#endif

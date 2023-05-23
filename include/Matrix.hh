#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
using std::vector;

//inherit from vector of doubles?
class Matrix{
	public:
		//set dimensionality
		Matrix();
		Matrix(int dim);
		Matrix(int x_dim, int y_dim);
		//constructor from inputs
		virtual ~Matrix();

		void InitRandom(double min, double max);	
		void SetDims(int x_dim, int y_dim);
		void SetEntry(double val, int i, int j);
		double GetEntry(int i, int j);
		double det();
		Matrix transpose();
		vector<int> GetDims();
		Matrix mult(const Matrix& mat);
		//update entries
		//add entries
		//clear
	private:
		int m_Xdim;
		int m_Ydim;
		vector<vector<double> m_entries;
		void Matrix::subMatrix(int mat[N][N], int temp[N][N], int p, int q, int n);

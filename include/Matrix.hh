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
		double at(int i, int j) const;
		void GetCofactor(Matrix& temp, int p, int q, int n);
		vector<int> GetDims() const;
		Matrix mult(const Matrix& mat) const;
		double det() const;
		double det(int n);
		bool square() const { return _square; }
		void transpose(const Matrix& mat);
		void adjoint(Matrix& adj);
		void invert(Matrix& inv);
		//clear
		void clear();
		//update entries
		//add entries
	private:
		int m_Xdim;
		int m_Ydim;
		vector<vector<double>> m_entries;
		//void subMatrix(int mat[N][N], int temp[N][N], int p, int q, int n);
		double m_det;
		bool _square;
};


#endif

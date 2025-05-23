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
		Matrix(double pt);
		Matrix(BayesPoint pt);
		Matrix(PointCollection pts);
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
		int Dim(int d) const{if(d == 0) return m_row; else if(d == 1) return m_col; else return -999;}
		//multiplies mat1*mat2 and stores result in this
		void mult(const Matrix& mat1, const Matrix& mat2);
		//multiplies mat by factor and stores in this
		void mult(const Matrix& mat, double factor);
		//add two matrices and store result in this
		void add(const Matrix& mat1, const Matrix& mat2);
		//add mat to this
		void add(const Matrix& mat);
		//subtract two matrices and store result in this
		void minus(const Matrix& mat1, const Matrix& mat2);
		//minus mat from this
		void minus(const Matrix& mat);
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
		PointCollection MatToPoints(vector<double> weights = {});
		void PointsToMat(PointCollection& pc);
		void PointToMat(const BayesPoint& pc);
		void mean(const PointCollection& data);
		void scatter(const PointCollection& data);

		bool empty() const{ if(m_row == 0 && m_col == 0) return true; else return false; }
		//clear
		void clear();
		void Print() const;
		void reset(){ clear(); InitEmpty(); }
		//sets this to be a shift matrix
		void PointToShift(const BayesPoint& pt){
			if(pt.Dim() != m_row) return;
			for(int i = 0; i < m_row; i++)
				for(int j = 0; j < m_col; j++)
					SetEntry(pt.at(i),i,j);
		}
		//sets this to be a scale matrix
		void PointToScale(const BayesPoint& pt){
			if(pt.Dim() != m_row) return;
			for(int i = 0; i < m_row; i++)
				SetEntry(pt.at(i),i,i);
		}
		//Frobenius inner product bw this and mat
		double FrobProd(const Matrix& mat){
			//check dims are the same
			if(m_row != mat.m_row || m_col != mat.m_col) return -999;
			double ret = 0;
			for(int i = 0 ; i < m_row; i++){
				for(int j = 0; j < m_col; j++){
					ret += m_entries[i][j]*mat.m_entries[i][j]; //for real-valued matrices
				}
			}
			return ret;
		}
	private:
		int m_row;
		int m_col;
		vector<vector<double>> m_entries;
		bool _square;
};


#endif

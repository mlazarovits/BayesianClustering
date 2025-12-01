#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include "PointCollection.hh"
#include <Eigen/Dense>
using std::vector;

class Matrix{
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW //to ensure proper alignment for heap allocation and STL containers
		//set dimensionality
		Matrix();
		//x_dim = dim (rows), y_dim = # points (cols)
		Matrix(int row, int col);
		Matrix(vector<double>& in);
		Matrix(double pt);
		Matrix(const BayesPoint& pt);
		Matrix(const PointCollection& pts);
		//rule of three
		//copy constructor
		Matrix(const Matrix& mat) = default;
		//copy assignment operator
		Matrix& operator=(const Matrix& mat) = default;
		//destructor
		virtual ~Matrix() = default;

		void InitRandom(double min = 0, double max = 1., unsigned long long seed = 123);
		//void InitRandomSym(double min = 0, double max = 1., unsigned long long seed = 123);
		//void InitRandomSymPosDef(double min = 0, double max = 1., unsigned long long seed = 333);
		void InitEmpty();
		//initializes identity matrix
		void InitIdentity();
		//fills this matrix with data pts from n-dim gaus
		void SampleNDimGaussian(const Matrix& mean, const Matrix& sigma, int Nsample);
		void SetDims(int row, int col);
		void SetEntry(double val, int i, int j);
		double at(int i, int j) const;
		//void GetCofactor(Matrix& temp, int p, int q, int n) const;
		int nRows() const;
		int nCols() const;
		int Dim(int d) const{if(d == 0) return _mat.rows(); else if(d == 1) return _mat.cols(); else return -999;}
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
		double det() const;
		bool square() const { return _square; }
		bool symmetric(double tol = 1e-10) const;
		void transpose(const Matrix& mat);
		void invert(const Matrix& mat);
		double trace();
		void cholesky(Matrix& L) const;
		void eigenCalc(vector<double>& vals, vector<Matrix>& vecs) const;
		void MatToPoints(PointCollection& pc, const vector<double>& weights = {}, const vector<int>& skipdims = {}, const vector<int>& idxs = {}) const;
		void MatToPoint(BayesPoint& pt, double weight = -999, int skipdim = -999, int idx = -999) const;
		void PointsToMat(const PointCollection& pc);
		void PointToMat(const BayesPoint& pc);
		void mean(const PointCollection& data);
		void scatter(const PointCollection& data);

		bool empty() const{ return (_mat.size() == 0); } 
		//clear
		void clear();
		void Print() const;
		void reset(){ InitEmpty(); }

	private:
		Eigen::MatrixXd _mat;
		bool _square;
};


#endif

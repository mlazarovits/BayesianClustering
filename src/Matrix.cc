#include "Matrix.hh"
#include "RandomSample.hh"
#include "PointCollection.hh"
#include <iostream>

using std::cout;
using std::endl;


Matrix::Matrix(){
	_square = true;
	_mat = Eigen::MatrixXd();
}


Matrix::Matrix(int row, int col){
	_mat = Eigen::MatrixXd(row, col);
	InitEmpty();
	if(_mat.rows() == _mat.cols()) _square = true;
}

Matrix::Matrix(vector<double>& in){
	_mat = Eigen::Map<Eigen::MatrixXd>(in.data(), (int)in.size(), 1);
} 


Matrix::Matrix(double pt){
	_mat = Eigen::MatrixXd(1,1);
	_mat(0,0) = pt;

}

Matrix::Matrix(const BayesPoint& pt){
	vector<double> vals;
	pt.Value(vals);
	_mat = Eigen::Map<Eigen::MatrixXd>(vals.data(), pt.Dim(), 1);
}

Matrix::Matrix(const PointCollection& pts){
	int row = pts.Dim();
	int col = pts.GetNPoints();

	std::vector<double> buffer(row * col);
	std::vector<double> pt;
	for(int j = 0; j < col; j++){
		pts.at(j).Value(pt);
		std::copy(pt.begin(), pt.end(), buffer.begin() + j * row);
	}
	_mat = Eigen::Map<Eigen::MatrixXd>(buffer.data(), row, col);

}


/*
//copy constructor
Matrix::Matrix(const Matrix& mat){
	_square = mat._square;
	_mat = mat._mat;
}	


Matrix::~Matrix(){ }

Matrix& Matrix::operator=(const Matrix& mat) {
    if (this != &mat) {           // protect against self-assignment
        _square = mat._square;
        _mat = mat._mat;          // Eigen assignment is safe
    }
    return *this;
}
*/

//creates a random matrix
void Matrix::InitRandom(double min, double max, unsigned long long seed){
	if(_mat.size() > 0){
		cout << "InitRandom Need dimensions to init." << endl;
		return;
	}
	RandomSample rs(seed);
	rs.SetRange(min, max);

	_mat = _mat.unaryExpr([&](double){return rs.SampleFlat();});

}


//creates an empty matrix
void Matrix::InitEmpty(){
	if(_mat.rows() == 0 && _mat.cols() == 0){
		cout << "InitEmpty - Need dimensions to init. - mat size "  << _mat.size() << " # rows " << _mat.rows() << " # cols " << _mat.cols() << endl;
		return;
	}
	_mat.setZero();
}

//creates the identity matrix
void Matrix::InitIdentity(){
	if(_mat.rows() == 0 && _mat.cols() == 0){
		cout << "InitIdentity - Need dimensions to init. - mat size "  << _mat.size() << " # rows " << _mat.rows() << " # cols " << _mat.cols() << endl;
		return;
	}
	if(_mat.rows() != _mat.cols()){ cout << "Matrix::InitEmpty - Non-square matrix." << endl; return;}


	_mat.setIdentity();

	_square = true;

}

void Matrix::SetDims(int row, int col){
	_mat = Eigen::MatrixXd(row, col);
	InitEmpty();
	if(_mat.rows() == _mat.cols()) _square = true;
}


void Matrix::SetEntry(double val, int i, int j){
	_mat(i,j) = val;
}


//n is current dim of cofactor
double Matrix::det() const{
	//unstable for large matrices, uses formulas for 2x2 and 3x3 matrices
	return _mat.determinant();

}


//invert mat and store in this
void Matrix::invert(const Matrix& mat){
	if(!mat.square()){
		cout << "Error: non-square matrix, cannot calculate inverse." << _mat.rows() << " " << _mat.cols() << endl;
		return;
	}
	_mat = mat._mat.inverse();
}




int Matrix::nRows() const{ return _mat.rows(); }
int Matrix::nCols() const{ return _mat.cols(); }

double Matrix::at(int i, int j) const{
	if(i > _mat.rows()-1 || j > _mat.cols()-1){
		cout << "Error: accessing element (" << i << "," << j << ") for matrix of dimension (" << _mat.rows() << "," << _mat.cols() << ")" << endl; 
		return -999; }
	return _mat(i,j);
}


//multiply mat by factor and store in this
void Matrix::mult(const Matrix& mat, double factor){
	_mat = mat._mat * factor;
}

//multiply two matrices and store result in this matrix
void Matrix::mult(const Matrix& mat1, const Matrix& mat2){
	_mat = mat1._mat * mat2._mat;
}



void Matrix::transpose(const Matrix& mat){
	_mat = mat._mat.transpose();
}
			
//clear
void Matrix::clear(){
	_mat = Eigen::MatrixXd();
}



bool Matrix::symmetric(double tol) const{
	return (_mat - _mat.transpose()).cwiseAbs().maxCoeff() < tol;

}

void Matrix::Print() const{
	cout << _mat << endl;
}




void Matrix::add(const Matrix& mat1, const Matrix& mat2){
	_mat = mat1._mat + mat2._mat;
}



void Matrix::add(const Matrix& mat){
	_mat += mat._mat;	
}



void Matrix::minus(const Matrix& mat){
	_mat -= mat._mat;
}


void Matrix::minus(const Matrix& mat1, const Matrix& mat2){
	_mat = mat1._mat - mat2._mat;


}


void Matrix::MatToPoint(BayesPoint& pt, double weight, int skipdim, int idx) const{
	if(_mat.cols() > 1){
		cout << "Error: cannot turn matrix of dimensions " << _mat.cols() << " x " << _mat.rows() << " into single point. Try MatToPoints(PointCollection pc)." << endl;
		return;
	}
	vector<double> val(_mat.data(), _mat.data() + _mat.size());
	pt = BayesPoint(val);
	if(weight != -999)
		pt.SetWeight(weight);
	if(skipdim != -999)
		pt.SetSkipDim(skipdim);
	if(idx != -999)
		pt.SetUserIdx(idx);
}



void Matrix::MatToPoints(PointCollection& pc, const vector<double>& weights, const vector<int>& skipdims, const vector<int>& idxs) const{
	pc.Clear(); 
	for(int j = 0; j < _mat.cols(); j++){
		Eigen::VectorXd col = _mat.col(j);
		vector<double> val(col.data(), col.data() + col.size());
		BayesPoint pt = BayesPoint(val);
		if(weights.size() == _mat.cols())
			pt.SetWeight(weights[j]);
		if(skipdims.size() == _mat.cols())
			pt.SetSkipDim(skipdims[j]);
		if(idxs.size() == _mat.cols())
			pt.SetUserIdx(idxs[j]);
		pc += pt;
	}
}

/*
//from Numerical Recipes 3rd Ed., Ch. 2.9
void Matrix::cholesky(Matrix& L) const{
	//check if matrix is square, symmetric
	L = Matrix(_mat.rows(), _mat.cols());
	if(!_square){
		cout << "Error: matrix is not square." << endl;
		return;
	}
	if(!symmetric()){
		cout << "Error: matrix is not symmetric." << endl;
		return;
	}
	double sum;
	int k;
	for(int i = 0; i < _mat.rows(); i++){
		for(int j = i; j < _mat.cols(); j++){
			//copy mat into L
			L.SetEntry(m_entries[i][j],i,j);
			for(sum = m_entries[i][j], k = i - 1; k >= 0; k--){
				sum -= L.at(i,k)*L.at(j,k);  
			}
			if(i == j){
				if(sum <= 0.0){
					cout << "mat not posdef. " << sum << endl;
					Print();	
					return;
				}
				L.SetEntry(sqrt(sum),i,i);
			}
			else{
				L.SetEntry(sum/L.at(i,i),j,i);
			}
		}
	} 
	//set upper-triangular part to 0 - L should be lower-triangular
	for(int i = 0; i < _mat.rows(); i++)
		for(int j = 0; j < i; j++)
			L.SetEntry(0.,j,i);
	
}

//from Numerical Recipes 3rd Ed., Ch. 7.4
//return Nsample random d-dim points (xmin, xmax) according to n-dim Gaussian distribuion
void Matrix::SampleNDimGaussian(const Matrix& mean, const Matrix& sigma, int Nsample){
	//need cholesky decomposition: takes sigma = LL^T
	//returns Nsample normally distributed n-dim points (one point = x)
	//x = L*y + mu for vectors x, y, mu and matrix L
	int d = mean.nRows();
	SetDims(d,Nsample);
	InitEmpty();
//	//aggregate points
	Matrix L;
	sigma.cholesky(L);
	RandomSample rs;
	vector<double> y;
	for(int n = 0; n < Nsample; n++){
		//y_i ~ N(0,1) - takes d independent samples from N(0,1)
		y = rs.SampleGaussian(0.,1.,d);
		Matrix mat_y(y);
		Matrix Ly_mu = Matrix(d,1);
		//L*y	
		Ly_mu.mult(L,mat_y);	
		//L*y + mu
		Ly_mu.add(mean);
		//add point to mat
		for(int j = 0; j < d; j++)
			m_entries[j][n] = Ly_mu.at(j,0);
	}
	cout << "Sampled " << Nsample << " " << d << "-dimensional points from a multidim Gaussian via cholesky decomposition." << endl;
}
*/

//fill vals + vecs vectors
void Matrix::eigenCalc(vector<double>& vals, vector<Matrix>& vecs) const{
	vals.clear();
	vecs.clear();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(_mat);
	if( eigensolver.info() != Eigen::Success){
		cout << "eigenCalc Error: not able to calculate eigenvalues/vectors." << endl;
		return;
	}
 
	vals.assign(eigensolver.eigenvalues().begin(), eigensolver.eigenvalues().begin() + eigensolver.eigenvalues().size());

	//columns are eigenvectors
	for(int d = 0; d < _mat.rows(); d++){
		vecs.push_back(Matrix(_mat.rows(), 1));
		vecs[d]._mat = eigensolver.eigenvectors().row(d);
	}
}



double Matrix::trace(){
	return _mat.trace();
}


//stores pc as mat in this
void Matrix::PointsToMat(const PointCollection& pc){
	int col = pc.GetNPoints();
	int row = pc.Dim();
	_mat = Eigen::MatrixXd(row, col);
	InitEmpty();
	for(int j = 0; j < col; j++){
		vector<double> in;
		pc.at(j).Value(in);
		_mat.col(j) = Eigen::Map<Eigen::VectorXd>(in.data(), in.size());
	}
}
//stores pc as mat in this
void Matrix::PointToMat(const BayesPoint& pc){
	int col = 1;
	int row = pc.Dim();
	_mat = Eigen::MatrixXd(row, col);
	InitEmpty();

	vector<double> in;
	pc.Value(in);
	_mat = Eigen::Map<Eigen::MatrixXd>(in.data(), row, col);

}



void Matrix::mean(const PointCollection& data){
	PointsToMat(data);
	_mat = _mat.rowwise().mean();
}




void Matrix::scatter(const PointCollection& data){
	int n = data.GetNPoints();
	int d = data.Dim();

	Matrix m = Matrix(d, 1);
	Matrix mT = Matrix(1, d);
	m.mean(data);
	mT.transpose(m);
	
	Matrix pt = Matrix(d,1);
	Matrix ptT = Matrix(1,d);
	Matrix pt_ptT = Matrix(d,d);
	for(int i = 0; i < n; i++){
		pt.PointToMat(data.at(i));
		ptT.transpose(pt);
		pt_ptT.mult(pt,ptT);
		add((*this),pt_ptT);
	}	
	
	pt_ptT.InitEmpty();
	pt_ptT.mult(m,mT);
	pt_ptT.mult(pt_ptT,n);

	minus(pt_ptT);

}

#include "Matrix.hh"
#include "RandomSample.hh"
#include "PointCollection.hh"
#include <iostream>
#include <Eigen/Dense>

using std::cout;
using std::endl;


Matrix::Matrix(){
	_square = true;
	m_row = 0;
	m_col = 0;
}


Matrix::Matrix(int row, int col){
	m_row = row;
	m_col = col;
	if(m_row == m_col) _square = true;
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(0.);
		}
	}		
}

Matrix::Matrix(vector<double> in){
	m_row = (int)in.size();
	m_col = 1;
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(in[i]);
		}
	}
} 


Matrix::Matrix(double pt){
	m_row = 1;
	m_col = 1;
	m_entries.push_back({});
	m_entries[0].push_back(pt);

}

Matrix::Matrix(BayesPoint pt){
	m_row = pt.Dim();
	m_col = 1;

	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(pt.at(i));
		}
	}
}

Matrix::Matrix(PointCollection pts){
	m_row = pts.Dim();
	m_col = pts.GetNPoints();
	InitEmpty();
	
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m_entries[i][j] = pts.at(j).at(i);
		}
	}

}



//copy constructor
Matrix::Matrix(const Matrix& mat){
	m_row = mat.GetDims()[0];
	m_col = mat.GetDims()[1];

	_square = mat.square();

	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back( mat.at(i,j) );
		}
	}
}	


Matrix::~Matrix(){
//	for(int i = 0; i < m_entries.size(); i++){
//		m_entries[i].clear();
//	}

}
//creates a random matrix
void Matrix::InitRandom(double min, double max, unsigned long long seed){
	if(m_entries.size() < 0){
		cout << "Need dimensions to init." << endl;
		return;
	}
	RandomSample rs(seed);
	//TODO: set range by data
	rs.SetRange(min, max);
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m_entries[i][j] = rs.SampleFlat();
		}
	}		
}

//creates a symmetric random matrix
void Matrix::InitRandomSym(double min, double max, unsigned long long seed){
	if(m_entries.size() < 0){
		cout << "Need dimensions to init." << endl;
		return;
	}
	RandomSample rs(seed);
	//TODO: set range by data
	rs.SetRange(min, max);
	if(!_square){
		cout << "Error: cannot initiate a symmetric matrix because dimensions provided were not equal (needs to be a square matrix)" << endl;
		return;
	}
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			if(i <= j){
				m_entries[i][j] = rs.SampleFlat();
			}
			else if(i > j){
				m_entries[i][j] = m_entries[j][i];
			}
		}
	}	
}


//creates a random symmetric positive definite matrix
void Matrix::InitRandomSymPosDef(double min, double max, unsigned long long seed){
	if(m_entries.size() < 0){
		cout << "Need dimensions to init." << endl;
		return;
	}
	if(!_square){
		cout << "Error: cannot initiate a symmetric matrix because dimensions provided were not equal (needs to be a square matrix). Rows: " << m_row << " cols: " << m_col << endl;
		return;
	}
	Matrix A = Matrix(m_row, m_col);
	Matrix A_T;// = Matrix(m_col, m_row);
	A.InitRandom(min,max,seed);
	A_T.transpose(A);
	
//	for(int i = 0; i < m_row; i++){
//		for(int j = 0; j < m_col; j++){
//			//cout << "i: " << i << " j: " << j << " A: " << A.at(i,j) << endl;
//			A_T.SetEntry(A.at(i,j),j,i);
//			//cout << "i: " << i << " j: " << j << " A: " << A.at(i,j) << " A_T: " << A_T.at(i,j) << endl;
//		}
//	}
	//Matrix sym(m_row, m_col);
	//sym.mult(A_T,A);
	//cout << "sym" << endl;
	//sym.Print();
	this->mult(A_T,A);
//	double val;
//	for(int i = 0; i < m_row; i++){
//	cout << "m_entries row " << i << " size: " << m_entries[i].size() << endl;
//		for(int j = 0; j < m_col; j++){
//			cout << "i: " << i << " j: " << j << " m_entries: " << m_entries[i][j] << endl;
//			val = 0;
//			for(int k = 0; k < m_row; k++){
//				val += A_T.at(i,k)*A.at(k,j);
//			}
//			m_entries[i][j] = val;
//		}
//	}
}

//creates an empty matrix
void Matrix::InitEmpty(){
	if(m_entries.size() > 0){
	//	cout << "Matrix already initialized: " << m_entries.size() << endl;
		return;
	}
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(0.);
		}
	}		
}

//creates the identity matrix
void Matrix::InitIdentity(){
	if(m_entries.size() < 0){
		cout << "Need dimensions to init." << endl;
		return;
	}
	if(m_row != m_col){ cout << "Matrix::InitEmpty - Non-square matrix." << endl; return;}
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			if(i == j)	
				m_entries[i][j] = 1.;
			else
				m_entries[i][j] = 0.;
		}
	}	
	_square = true;

}

void Matrix::SetDims(int row, int col){
	m_row = row;
	m_col = col;
	if(m_row == m_col) _square = true;
	if(m_entries.size() < 1)
		InitEmpty();
}


void Matrix::SetEntry(double val, int i, int j){
	m_entries[i][j] = val;

}


void Matrix::GetCofactor(Matrix& temp, int p, int q, int n) const{
	int i = 0, j = 0;
	//int n = temp.GetDims()[0];
	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  Copying into temporary matrix only those elements which are not in given row and column
			if (row != p && col != q) {
				temp.SetEntry(m_entries[row][col],i,j++);
				// Row is filled, so increase row index and reset col index
				if (j == n - 1) {
				    j = 0;
				    i++;
				}
			}
		}
	}
}


//n is current dim of cofactor
double Matrix::det(int n) const{
	double D = 0; // Initialize result
	if(!_square){
		cout << "Error: non-square matrix, cannot calculate determinant." << m_row << " " << m_col << endl;
		return -999;
	}
	if( n == 0){
		n = m_row;
	}

	//  Base case : if matrix contains single element
	if (n == 1)
	    return m_entries[0][0];
	
	Matrix temp(n,n); // To store cofactors
	
	int sign = 1; // To store sign multiplier
	
	// Iterate for each element of first row
	for (int f = 0; f < n; f++) {
	    // Getting Cofactor of m_entries[0][f]
	    GetCofactor(temp, 0, f, n);
	    D += sign * m_entries[0][f] * temp.det(n - 1);
	
	    // terms are to be added with alternate sign
	    sign = -sign;
	}
	return D;

}

// Function to get adjoint of m_entries[N][N] in adj[N][N].
void Matrix::adjoint(const Matrix& mat){
	if(m_row == 0 && m_col == 0){
		SetDims(mat.GetDims()[0],mat.GetDims()[1]);
		InitEmpty();
	}
	if(m_row != mat.GetDims()[0] || m_col != mat.GetDims()[1]){
		cout << "Error: dimensions of this matrix do not match mat (matrix given for adjoint calculation." << endl;
		return;
	}
	if(!_square){
		cout << "Error: non-square matrix, cannot calculate adjoint." << m_row << " " << m_col << endl;
		return;
	}

	if (m_row == 1) {
	    m_entries[0][0] = 1;
	    return;
	}
	
	// temp is used to store cofactors of m_entries[][]
	int sign = 1;
	int N = m_row;
	Matrix temp = Matrix(N,N);
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			// Get cofactor of m_entries[i][j]
			mat.GetCofactor(temp, i, j, N);
			// sign of adj[j][i] positive if sum of row and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;
			
			// Interchanging rows and columns to get the transpose of the cofactor matrix
			m_entries[j][i] = (sign) * (temp.det(N - 1));
		}
	}
}



//invert mat and store in this
void Matrix::invert(const Matrix& mat){
	if(!mat.square()){
		cout << "Error: non-square matrix, cannot calculate inverse." << m_row << " " << m_col << endl;
		return;
	}
	vector<int> dims = mat.GetDims();
	SetDims(dims[0],dims[1]);	
	InitEmpty();
	double det = mat.det(dims[0]);
	if (det == 0) {
	    cout << "Singular matrix, can't find its inverse" << endl;
	    return;
	}
	Eigen::MatrixXd m(m_row, m_col);

	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m(i,j) = mat.at(i,j);
		}
	}

	Eigen::MatrixXd minv = m.inverse();

	for (int i = 0; i < dims[0]; i++)
	    for (int j = 0; j < dims[1]; j++)
	        m_entries[i][j] = minv(i,j);

}








vector<int> Matrix::GetDims() const{
	return {m_row, m_col};
}

double Matrix::at(int i, int j) const{
	if(i > m_row-1 || j > m_col-1){
		cout << "Error: accessing element (" << i << "," << j << ") for matrix of dimension (" << m_row << "," << m_col << ")" << endl; 
		return -999; }
	return m_entries[i][j];
}


//multiply mat by factor and store in this
void Matrix::mult(const Matrix& mat, double factor){
	//SetDims(mat.GetDims()[0], mat.GetDims()[1]);
	//InitEmpty();
	//in case mat = this (ie m = c*m)
	Matrix tmp(mat.GetDims()[0],mat.GetDims()[1]);
	for(int i = 0; i < m_row; i++)
		for(int j = 0; j < m_col; j++)
			tmp.SetEntry(factor*mat.at(i,j),i,j);
			//m_entries[i][j] = factor*mat.at(i,j);
	*this = tmp;
}

//multiply two matrices and store result in this matrix
void Matrix::mult(const Matrix& mat1, const Matrix& mat2){
	double val;
	vector<int> dims1 = mat1.GetDims();
	vector<int> dims2 = mat2.GetDims();
	//clear current matrix to be filled with mutliplied matrix
	if(dims1[1] != dims2[0]){
		cout << "Matrix::mult error: matrices have incompatible dimensions. mat1: " << dims1[0] << " x " << dims1[1] << " mat2: " << dims2[0] << " x " << dims2[1] << endl;
		return;
	}
	if(m_row != 0 and m_col != 0){
		if(m_row != dims1[0] or m_col != dims2[1]){
			cout << "Error: this matrix must be " << dims1[0] << " x " << dims2[1] << " dims (mult() doesn't change dimensions of this matrix)." << endl;
			return;
		}
	}
	//else{
	//	//SetDims(dims1[0],dims2[1]);
	//	//InitEmpty();
	//}
	//need to create temp matrix in case one of the passed matrices is this (ie m = A*m)
	Matrix tmp(dims1[0],dims2[1]);
		
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims2[1]; j++){
			val = 0;
			for(int k = 0; k < dims2[0]; k++){
				val += mat1.at(i,k) * mat2.at(k,j);
			}
			//m_entries[i][j] = val;
			tmp.SetEntry(val,i,j);
		}
	}
	*this = tmp;
}



void Matrix::transpose(const Matrix& mat){
	clear();
	SetDims(mat.GetDims()[1], mat.GetDims()[0]);
	InitEmpty();
	
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m_entries[i][j] = mat.at(j,i);
		}
	}
}
			
//clear
void Matrix::clear(){
	for(int i = 0; i < m_row; i++){
		m_entries[i].clear();
	}
	m_entries.clear();
}



bool Matrix::symmetric() const{
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			if(i < j) break; //only need to look at one half of the matrix
			if(m_entries[i][j] != m_entries[j][i]){
				return false;
			}
		}
	}
	return true;

}

void Matrix::Print() const{
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			cout << m_entries[i][j] << " ";
		}
		cout << endl;
	}
}




void Matrix::add(const Matrix& mat1, const Matrix& mat2){
	//check dims are compatible
	vector<int> dims1 = mat1.GetDims();
	vector<int> dims2 = mat2.GetDims();
	if(dims1[0] != dims2[0] || dims1[1] != dims2[1]){
		cout << "Error: matrix dimensions are not compatible for addition: " << dims1[0] << "x" << dims1[1] << " and " << dims2[0] << "x" << dims2[1] << endl;
		return;
	}
	SetDims(dims1[0],dims1[1]);
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims1[1]; j++){
			m_entries[i][j] = mat1.at(i,j) + mat2.at(i,j);
		}
	}
}



void Matrix::add(const Matrix& mat){
	//check dims are compatible
	//cout << "MATRIX::ADD - START" << endl;
	vector<int> dims1 = mat.GetDims();
	if(dims1[0] != m_row || dims1[1] != m_col){
		cout << "Error: matrix dimensions are not compatible for addition: " << dims1[0] << "x" << dims1[1] << " and " << m_row << "x" << m_col << endl;
		return;
	}
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims1[1]; j++){
			//cout << "i " << i << " j " << j << " og entry " << m_entries[i][j] << " + " << mat.at(i,j) << endl;
			m_entries[i][j] += mat.at(i,j);
			//cout << "i " << i << " j " << j << " new entry " << m_entries[i][j] << endl;
		}
	}
	//cout << "MATRIX::ADD - END" << endl;
	
}



void Matrix::minus(const Matrix& mat){
	Matrix sub;
	sub.mult(mat,-1);
	add(sub);
}


void Matrix::minus(const Matrix& mat1, const Matrix& mat2){
	Matrix sub2;
	sub2.mult(mat2, -1);
	add(mat1, sub2);



}


PointCollection Matrix::MatToPoints(vector<double> weights){
	PointCollection pc;
	for(int j = 0; j < m_col; j++){
		vector<double> val;
		for(int i = 0; i < m_row; i++){
			val.push_back(m_entries[i][j]);
		}
		BayesPoint pt = BayesPoint(val);
		if(weights.size() == m_col)
			pt.SetWeight(weights[j]);
		pc += pt;
	}
	return pc;
}

//from Numerical Recipes 3rd Ed., Ch. 2.9
Matrix Matrix::cholesky(){
	//check if matrix is square, symmetric
	Matrix L = Matrix(m_row, m_col);
	if(!_square){
		cout << "Error: matrix is not square." << endl;
		return L;
	}
	if(!symmetric()){
		cout << "Error: matrix is not symmetric." << endl;
		return L;
	}
	double sum;
	int k;
	for(int i = 0; i < m_row; i++){
		for(int j = i; j < m_col; j++){
			//copy mat into L
			L.SetEntry(m_entries[i][j],i,j);
			for(sum = m_entries[i][j], k = i - 1; k >= 0; k--){
				sum -= L.at(i,k)*L.at(j,k);  
			}
			if(i == j){
				if(sum <= 0.0){
					cout << "mat not posdef. " << sum << endl;
					Print();	
					return L;
				}
				L.SetEntry(sqrt(sum),i,i);
			}
			else{
				L.SetEntry(sum/L.at(i,i),j,i);
			}
		}
	} 
	//set upper-triangular part to 0 - L should be lower-triangular
	for(int i = 0; i < m_row; i++)
		for(int j = 0; j < i; j++)
			L.SetEntry(0.,j,i);
	return L;
	
}

//from Numerical Recipes 3rd Ed., Ch. 7.4
//return Nsample random d-dim points (xmin, xmax) according to n-dim Gaussian distribuion
void Matrix::SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample){
	//need cholesky decomposition: takes sigma = LL^T
	//returns Nsample normally distributed n-dim points (one point = x)
	//x = L*y + mu for vectors x, y, mu and matrix L
	int d = mean.GetDims()[0];
	SetDims(d,Nsample);
	InitEmpty();
//	//aggregate points
	Matrix L = sigma.cholesky();	
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


//fill vals + vecs vectors
void Matrix::eigenCalc(vector<double>& vals, vector<Matrix>& vecs){
	vals.clear();
	vecs.clear();

	Eigen::MatrixXd m(m_row, m_col);

	for(int i = 0; i < m_row; i++){
		vecs.push_back( Matrix(m_row,1) );
		for(int j = 0; j < m_col; j++){
			m(i,j) = m_entries[i][j];
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m);
	if( eigensolver.info() != Eigen::Success){
		cout << "eigenCalc Error: not able to calculate eigenvalues/vectors." << endl;
		return;
	}
 
	//columns are eigenvectors
	for(int d = 0; d < m_row; d++){
		for(int j = 0; j < m_col; j++){
			vecs[d].SetEntry(eigensolver.eigenvectors()(d,j),j,0);
		}
		vals.push_back(eigensolver.eigenvalues()(d));
	}

}



double Matrix::trace(){
	double tr = 0;
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			if(i == j)
				tr += m_entries[i][j];
		}
	}
	return tr;
}


//stores pc as mat in this
void Matrix::PointsToMat(PointCollection& pc){
	m_col = pc.GetNPoints();
	m_row = pc.Dim();

	InitEmpty();
	
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m_entries[i][j] = pc.at(j).at(i);
		}
	}
}
//stores pc as mat in this
void Matrix::PointToMat(const BayesPoint& pc){
	m_col = 1;
	m_row = pc.Dim();

	InitEmpty();
	
	for(int i = 0; i < m_row; i++){
			m_entries[i][0] = pc.at(i);
	}
}



void Matrix::mean(const PointCollection& data){
	int n = data.GetNPoints();
	int d = data.Dim();
	m_col = 1;
	m_row = d;
	InitEmpty();
	for(int j = 0; j < d; j++){
		for(int i = 0; i < n; i++){
			m_entries[j][0] += data.at(i).at(j); 
		}
	}
	for(int j = 0; j < d; j++) m_entries[j][0] = m_entries[j][0]/double(n);
}




void Matrix::scatter(const PointCollection& data){
	int n = data.GetNPoints();
	int d = data.Dim();
	m_row = d;
	m_col = d;

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

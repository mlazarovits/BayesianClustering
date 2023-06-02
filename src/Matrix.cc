#include "Matrix.hh"
#include "RandomSample.hh"
#include "PointCollection.hh"
#include <iostream>

using std::cout;
using std::endl;


Matrix::Matrix(){
	_square = true;
	_id = false;
	m_row = 0;
	m_col = 0;
}


Matrix::Matrix(int row, int col){
	m_row = row;
	m_col = col;
	_id = false;
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
	_id = false;
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(in[i]);
		}
	}
} 


//constructor from inputs

Matrix::~Matrix(){
//	for(int i = 0; i < m_entries.size(); i++){
//		m_entries[i].clear();
//	}

}

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
void Matrix::InitEmpty(){
	if(m_entries.size() > 0){
	//	cout << "Matrix already initialized." << endl;
		return;
	}
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(0.);
		}
	}		

}


void Matrix::InitIdentity(){
	if(m_entries.size() < 0){
		cout << "Need dimensions to init." << endl;
		return;
	}
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			if(i == j)	
				m_entries[i][j] = 1.;
			else
				m_entries[i][j] = 0.;
		}
	}	
	_id = true;

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
	//det(I) = 1
	if(_id){
		return 1.;
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


void Matrix::invert(const Matrix& mat){
	// Find determinant of m_entries[][]
	if(!mat.square()){
		cout << "Error: non-square matrix, cannot calculate inverse." << m_row << " " << m_col << endl;
		return;
	}
	vector<int> dims = mat.GetDims();
	SetDims(dims[0],dims[1]);	
	InitEmpty();
	double det = mat.det(dims[0]);
	//inverse of identity is itself
	if(mat.identity()){
		InitIdentity();
		return;
	}
	if (det == 0) {
	    cout << "Singular matrix, can't find its inverse" << endl;
	    return;
	}
	
	// Find adjoint
	Matrix adj = Matrix(dims[0], dims[1]);
	adj.adjoint(mat);
//	adj.Print();	
	// Find Inverse using formula "inverse(A) =
	// adj(A)/det(A)"
	for (int i = 0; i < dims[0]; i++)
	    for (int j = 0; j < dims[1]; j++)
	        m_entries[i][j] = adj.at(i,j)/det;
}




/*
double Matrix::det(){
	if(!_square){
		cout << "Matrix not square, determinant cannot be found." << m_row << " " << m_col << endl;
		return 0.;
	}
	int det = 0;
	Matrix submatrix = Matrix(m_row, m_col);
	if (n == 2)
	return ((m_entries[0][0] * m_entries[1][1]) - (m_entries[1][0] * m_entries[0][1]));
	else {
	   for (int x = 0; x < n; x++) {
	      int subi = 0;
	      for (int i = 1; i < n; i++) {
	         int subj = 0;
	         for (int j = 0; j < n; j++) {
	            if (j == x)
	            continue;
	            submatrix[subi][subj] = m_entries[i][j];
	            subj++;
	         }
	         subi++;
	      }
	      det = det + (pow(-1, x) * matrix[0][x] * determinant( submatrix, n - 1 ));
	   }
	}
	return det;
}
*/

vector<int> Matrix::GetDims() const{
	return {m_row, m_col};
}

double Matrix::at(int i, int j) const{
	return m_entries[i][j];
}


//multiply mat by factor and store in this
void Matrix::mult(const Matrix& mat, double factor){
	SetDims(mat.GetDims()[0], mat.GetDims()[1]);
	InitEmpty();
	for(int i = 0; i < m_row; i++)
		for(int j = 0; j < m_col; j++)
			m_entries[i][j] = factor*mat.at(i,j);
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
			cout << "Error: this matrix must be " << dims1[0] << " x " << dims2[1] << " dims." << endl;
			return;
		}
	}
	else{
		SetDims(dims1[0],dims2[1]);
		InitEmpty();
	}
	
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims2[1]; j++){
			val = 0;
			for(int k = 0; k < dims2[0]; k++){
				val += mat1.at(i,k) * mat2.at(k,j);
			}
			m_entries[i][j] = val;
		}
	}
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
		cout << "Error: matrix dimensions are not compatible for addition." << endl;
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
	vector<int> dims1 = mat.GetDims();
	if(dims1[0] != m_row || dims1[1] != m_col){
		cout << "Error: matrix dimensions are not compatible for addition." << endl;
		return;
	}
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims1[1]; j++){
			m_entries[i][j] += mat.at(i,j);
		}
	}
	
}


PointCollection Matrix::MatToPoints(){
	PointCollection pc;
	for(int j = 0; j < m_col; j++){
		vector<double> val;
		for(int i = 0; i < m_row; i++){
			val.push_back(m_entries[i][j]);
		}
		Point pt = Point(val);
		pc += pt;
	}
	return pc;
}

//from Numerical Recipes 3rd Ed., Ch. 2.9
Matrix Matrix::cholesky(){
	//check if matrix is square, posdef + symmetric
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
	rs.SetRange(0.,1.);
	vector<double> y;
	for(int n = 0; n < Nsample; n++){
		//y_i ~ N(0,1)
		y = rs.SampleGaussian(0.,1.,d);
		Matrix mat_y(y);
		//L*y	
		Matrix Ly_mu = Matrix(d,1);
		Ly_mu.mult(L,mat_y);	
		//L*y + mu
		Ly_mu.add(mean);
		//add point to mat
		for(int j = 0; j < d; j++)
			m_entries[j][n] = Ly_mu.at(j,0);
	}
	cout << "Sampled " << Nsample << " " << d << "-dimensional points from a multidim Gaussian via cholesky decomposition." << endl;
}




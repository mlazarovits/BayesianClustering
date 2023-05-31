#include "Matrix.hh"
#include "RandomSample.hh"
#include <iostream>

using std::cout;
using std::endl;


Matrix::Matrix(){
	m_det = -999;
	_square = false;
}

Matrix::Matrix(int dim){
	m_det = -999;
	_square = true;
	m_Xdim = dim;
	m_Ydim = dim;
	for(int i = 0; i < m_Xdim; i++){
		vector<double> vec;
		m_entries.push_back(vec);
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i].push_back(0.);
		}
	}		
}

Matrix::Matrix(int x_dim, int y_dim){
	m_det = -999;
	m_Xdim = x_dim;
	m_Ydim = y_dim;
	if(m_Xdim == m_Ydim) _square = true;
	for(int i = 0; i < m_Xdim; i++){
		vector<double> vec;
		m_entries.push_back(vec);
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i].push_back(0.);
		}
	}		
}

Matrix::Matrix(vector<double> in){
	m_Xdim = (int)in.size();
	m_Ydim = 1;

	m_entries.push_back({});
	for(int i = 0; i < m_Xdim; i++)
		m_entries[0].push_back(in[i]);
	
} 


//constructor from inputs

Matrix::~Matrix(){
	for(int i = 0; i < m_Xdim; i++){
		m_entries[i].clear();
	}
}

void Matrix::InitRandom(double min, double max){
	RandomSample rs;
	//TODO: set range by data
	rs.SetRange(min, max);
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i][j] = rs.SampleFlat();
		}
	}		
}

void Matrix::InitRandomSym(double min, double max){
	RandomSample rs;
	//TODO: set range by data
	rs.SetRange(min, max);
	if(!_square){
		cout << "Error: cannot initiate a symmetric matrix because dimensions provided were not equal (needs to be a square matrix)" << endl;
		return;
	}
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			if(i <= j){
				m_entries[i][j] = rs.SampleFlat();
			}
			else if(i > j){
				m_entries[i][j] = m_entries[j][i];
			}
		}
	}	
}
void Matrix::InitRandomSymPosDef(double min, double max){
	if(!_square){
		cout << "Error: cannot initiate a symmetric matrix because dimensions provided were not equal (needs to be a square matrix)" << endl;
		return;
	}
	Matrix A = Matrix(m_Xdim, m_Ydim);
	Matrix A_T = Matrix(m_Xdim, m_Ydim);
	A.InitRandom(min,max);
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			A_T.SetEntry(A.at(i,j),j,i);
		}
	}
	Matrix sym(m_Xdim, m_Ydim);
	sym.mult(A_T,A);
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i][j] = sym.at(i,j);
		}
	}
}
void Matrix::SetDims(int x_dim, int y_dim){
	m_Xdim = x_dim;
	m_Ydim = y_dim;
	if(m_Xdim == m_Ydim) _square = true;
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
	//if det has already been calculated, return result
//	if(m_det != -999){
//		return m_det;
//	}
	double D = 0; // Initialize result
	if(!_square){
		cout << "Error: non-square matrix, cannot calculate determinant." << m_Xdim << " " << m_Ydim << endl;
		return -999;
	}
	if( n == 0){
		n = m_Xdim;
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
//	m_det = D;	
	return D;

}

// Function to get adjoint of m_entries[N][N] in adj[N][N].
void Matrix::adjoint(Matrix& adj){
	if(!_square){
		cout << "Error: non-square matrix, cannot calculate adjoint." << m_Xdim << " " << m_Ydim << endl;
		return;
	}

	if (m_Xdim == 1) {
	    adj.SetEntry(1,0,0);
	    return;
	}
	
	// temp is used to store cofactors of m_entries[][]
	int sign = 1;
	int N = m_Xdim;
	Matrix temp = Matrix(N,N);
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			// Get cofactor of m_entries[i][j]
			GetCofactor(temp, i, j, N);
			
			// sign of adj[j][i] positive if sum of row and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;
			
			// Interchanging rows and columns to get the transpose of the cofactor matrix
			adj.SetEntry((sign) * (temp.det(N - 1)), j, i);
		}
	}
}


void Matrix::invert(Matrix& mat){
	// Find determinant of m_entries[][]
	if(!mat.square()){
		cout << "Error: non-square matrix, cannot calculate inverse." << m_Xdim << " " << m_Ydim << endl;
		return;
	}
	vector<int> dims = mat.GetDims();
	this->SetDims(dims[0],dims[1]);
	int det = mat.det(dims[0]);
	if (det == 0) {
	    cout << "Singular matrix, can't find its inverse" << endl;
	    return;
	}
	
	// Find adjoint
	Matrix adj = Matrix(dims[0], dims[1]);
	adjoint(adj);
	
	// Find Inverse using formula "inverse(A) =
	// adj(A)/det(A)"
	for (int i = 0; i < dims[0]; i++)
	    for (int j = 0; j < dims[1]; j++)
	        this->SetEntry(adj.at(i,j) / float(det), i, j);
	
}




/*
double Matrix::det(){
	if(!_square){
		cout << "Matrix not square, determinant cannot be found." << m_Xdim << " " << m_Ydim << endl;
		return 0.;
	}
	int det = 0;
	Matrix submatrix = Matrix(m_Xdim, m_Ydim);
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
	return {m_Xdim, m_Ydim};
}

double Matrix::at(int i, int j) const{
	return m_entries[i][j];
}


void Matrix::mult(const Matrix& mat1, const Matrix& mat2){
	double val;
	vector<int> dims1 = mat1.GetDims();
	vector<int> dims2 = mat2.GetDims();
	//clear current matrix to be filled with mutliplied matrix
	if(dims1[1] != dims2[0]){
		cout << "Error: matrices have incompatible dimensions. " << dims1[1] << " " << dims2[0] << endl;
		return;
	}
	clear();
	SetDims(dims1[0],dims2[1]);
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
	this->clear();
	this->SetDims(m_Ydim, m_Xdim);

	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			this->SetEntry(mat.at(j,i),i,j);
		}
	}
}
			
//clear
void Matrix::clear(){
	for(int i = 0; i < m_Xdim; i++){
		m_entries[i].clear();
	}
	m_entries.clear();
}

Matrix Matrix::cholesky(){
	//check if matrix is square, posdef + symmetric
	Matrix L = Matrix(m_Xdim, m_Ydim);
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
	for(int i = 0; i < m_Xdim; i++){
	//	cout << "i: " << i << endl;
		for(int j = i; j < m_Ydim; j++){
			//copy mat into L
			L.SetEntry(m_entries[i][j],i,j);
		//	cout << " 	j: " << j << endl;
			for(sum = m_entries[i][j], k = i - 1; k >= 0; k--){
		//cout << "		k: " << k << endl;
		//cout << "i: " << i << " j: " << j << " k: " << k << endl;
				sum -= L.at(i,k)*L.at(j,k);  
			}
			if(i == j){
				if(sum <= 0.0){
					cout << "mat not posdef. " << sum << endl;
					print();	
					return L;
				}
				L.SetEntry(sqrt(sum),i,i);
			}
			else{
				L.SetEntry(sum/L.at(i,i),j,i);
			}
			//	cout << "	sum: " << sum << endl;
		//cout << "L_ji = " << m_entries[j][i] << " for j: " << j << " i: " << i << endl;
		}
	} 
	//set upper-triangular part to 0 - L should be lower-triangular
	for(int i = 0; i < m_Xdim; i++)
		for(int j = 0; j < i; j++)
			L.SetEntry(0.,j,i);
	//print();
	return L;
	
}


bool Matrix::symmetric() const{
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			if(i < j) break; //only need to look at one half of the matrix
			if(m_entries[i][j] != m_entries[j][i]){
				return false;
			}
		}
	}
	return true;

}




void Matrix::print() const{
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
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
	if(dims1[0] != m_Xdim || dims1[1] != m_Ydim){
		cout << "Error: matrix dimensions are not compatible for addition." << endl;
		return;
	}
	for(int i = 0; i < dims1[0]; i++){
		for(int j = 0; j < dims1[1]; j++){
			m_entries[i][j] += mat.at(i,j);
		}
	}
	
}


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

//constructor from inputs

Matrix::~Matrix(){
	for(int i = 0; i < m_Xdim; i++){
		m_entries[i].clear();
	}
}

void Matrix::InitRandom(double min = 0., double max = 1.){
	RandomSample rs;
	//TODO: set range by data
	rs.SetRange(min, max);
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i][j] = rs.SampleFlat();
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


void Matrix::GetCofactor(Matrix& temp, int p, int q, int n){
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
double Matrix::det(int n = 0){
	//if det has already been calculated, return result
	if(m_det != -999){
		return m_det;
	}
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
	m_det = D;	
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


Matrix Matrix::mult(const Matrix& mat) const{
	double val;
	vector<int> dims = mat.GetDims();
	Matrix ret(m_Xdim, dims[1]);
	if(dims[0] != m_Ydim){
		cout << "Error: matrices have incompatible dimensions. " << dims[0] << " " << m_Ydim << endl;
		return ret;
	}
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			val = 0;
			for(int k = 0; k < dims[0]; k++){
				val += m_entries[i][k] * mat.at(k,j);
			}
			ret.SetEntry(val,i,j);
		}
	}
	return ret;
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
			
void Matrix::clear(){
	for(int i = 0; i < m_Xdim; i++){
		m_entries[i].clear();
	}
	m_entries.clear();
}




//update entries
//add entries
//clear


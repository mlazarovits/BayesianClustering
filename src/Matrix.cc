#include "Matrix.hh"
#include "RandomSample.hh"
#include "PointCollection.hh"
#include <iostream>

using std::cout;
using std::endl;


Matrix::Matrix(){
cout << "empty ctr" << endl;
	_square = true;
	m_row = 0;
	m_col = 0;
cout << "empty ctr - end" << endl;
}


Matrix::Matrix(int row, int col){
cout << "matrix ctr" << endl;
	m_row = row;
	m_col = col;
	if(m_row == m_col) _square = true;
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(0.);
			cout << "	i: " << i << " j: " << j << endl;
		}
	}		
cout << "matrix ctr - end" << endl;
}

Matrix::Matrix(vector<double> in){
cout << "matrix ctr" << endl;
	m_row = 1;
	m_col = (int)in.size();

	m_entries.push_back({});
	for(int j = 0; j < m_col; j++){
		m_entries[0].push_back(in[j]);
			cout << "	j: " << j << endl;
	}
cout << "matrix ctr - end" << endl;
} 


//constructor from inputs

Matrix::~Matrix(){
//	for(int i = 0; i < m_entries.size(); i++){
//		m_entries[i].clear();
//	}

}

void Matrix::InitRandom(double min, double max){
	cout << "InitRandom" << endl;
	RandomSample rs;
	//TODO: set range by data
	rs.SetRange(min, max);
	for(int i = 0; i < m_row; i++){
		for(int j = 0; j < m_col; j++){
			m_entries[i][j] = rs.SampleFlat();
			cout << "	i: " << i << " j: " << j << endl;
		}
	}		
	cout << "InitRandom - end" << endl;
}

void Matrix::InitRandomSym(double min, double max){
	RandomSample rs;
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
void Matrix::InitRandomSymPosDef(double min, double max){
	cout << "InitRandomSymPosDef" << endl;
	if(!_square){
		cout << "Error: cannot initiate a symmetric matrix because dimensions provided were not equal (needs to be a square matrix). Rows: " << m_row << " cols: " << m_col << endl;
		return;
	}
cout << "m_entries size: " << m_entries.size() << endl;
	Matrix A = Matrix(m_row, m_col);
	Matrix A_T;
	A.InitRandom(min,max);
	A_T.transpose(A);
//	for(int i = 0; i < m_row; i++){
//		for(int j = 0; j < m_col; j++){
//			//cout << "i: " << i << " j: " << j << " A: " << A.at(i,j) << endl;
//			A_T.SetEntry(A.at(i,j),j,i);
//			//cout << "i: " << i << " j: " << j << " A: " << A.at(i,j) << " A_T: " << A_T.at(i,j) << endl;
//		}
//	}
	cout << "A" << endl;
	//A.print();
	cout << "A_T" << endl;
//	A_T.print();
	//Matrix sym(m_row, m_col);
	//sym.mult(A_T,A);
	//cout << "sym" << endl;
	//sym.print();
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
	cout << "InitRandomSymPosDef - end" << endl;
}
void Matrix::InitEmpty(){
	if(m_entries.size() > 0){
		cout << "Matrix already initialized." << endl;
		return;
	}
	for(int i = 0; i < m_row; i++){
		m_entries.push_back({});
		for(int j = 0; j < m_col; j++){
			m_entries[i].push_back(0.);
		}
	}		

}


void Matrix::SetDims(int row, int col){
	m_row = row;
	m_col = col;
	if(m_row == m_col) _square = true;
	//init empty matrix
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
void Matrix::adjoint(Matrix& adj){
	if(!_square){
		cout << "Error: non-square matrix, cannot calculate adjoint." << m_row << " " << m_col << endl;
		return;
	}

	if (m_row == 1) {
	    adj.SetEntry(1,0,0);
	    return;
	}
	
	// temp is used to store cofactors of m_entries[][]
	int sign = 1;
	int N = m_row;
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
		cout << "Error: non-square matrix, cannot calculate inverse." << m_row << " " << m_col << endl;
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


void Matrix::mult(const Matrix& mat1, const Matrix& mat2){
	double val;
	vector<int> dims1 = mat1.GetDims();
	vector<int> dims2 = mat2.GetDims();
	//clear current matrix to be filled with mutliplied matrix
	if(dims1[1] != dims2[0]){
		cout << "Error: matrices have incompatible dimensions. " << dims1[1] << " " << dims2[0] << endl;
		return;
	}
	if(m_row != 0 and m_col != 0){
		if(m_row != dims1[0] or m_col != dims2[1]){
			cout << "Error: this matrix must be " << dims1[0] << " x " << dims2[1] << " dims." << endl;
			return;
		}
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
	SetDims(m_col, m_row);
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

void Matrix::print() const{
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
	for(int i = 0; i < m_row; i++){
		Point pt = Point(m_entries[i]);
		pc += pt;
		
	}
	return pc;
}

Matrix Matrix::cholesky(){
cout << "cholesky" << endl;
	cout << "rows: " << m_row << " cols: " << m_col << endl;
	//check if matrix is square, posdef + symmetric
	Matrix L = Matrix(m_row, m_col);
	return L;
	cout << "made L" << endl;
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
	cout << "start for loops" << endl;
	for(int i = 0; i < m_row; i++){
		for(int j = i; j < m_col; j++){
			cout << "1 - i: " << i << " j: " << j << " m_entries: " << m_entries[i][j] << endl;
			//copy mat into L
			L.SetEntry(m_entries[i][j],i,j);
			for(sum = m_entries[i][j], k = i - 1; k >= 0; k--){
				sum -= L.at(i,k)*L.at(j,k);  
			}
			if(i == j){
				if(sum <= 0.0){
					cout << "mat not posdef. " << sum << endl;
					print();	
					return L;
				}
				cout << "i: " << i << " j: " << j << " L: " << sqrt(sum) << endl;
				L.SetEntry(sqrt(sum),i,i);
			}
			else{
				cout << "i: " << i << " j: " << j << " L: " << sqrt(sum)/L.at(i,i) << endl;
				L.SetEntry(sum/L.at(i,i),j,i);
			}
		}
	} 
	//set upper-triangular part to 0 - L should be lower-triangular
	for(int i = 0; i < m_row; i++)
		for(int j = 0; j < i; j++)
			L.SetEntry(0.,j,i);
cout << "cholesky - end" << endl;
	return L;
	
}

//return Nsample random d-dim points (xmin, xmax) according to n-dim Gaussian distribuion
void Matrix::SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample){
	//need cholesky decomposition: takes sigma = LL^T
	//returns Nsample normally distributed n-dim points (one point = x)
	//x = L*y + mu for vectors x, y, mu and matrix L
	cout << "SampleNDimGaussian" << endl;
	int d = mean.GetDims()[0];
	//cout << "dim: " << d << endl;
//	
	//SetDims(d,Nsample);
//	cout << "m_row: " << m_row << endl;
	//InitEmpty();
//cout << "set dims" << endl;
//	//aggregate points
//	cout << "Sampling " << Nsample << " " << d << "-dim data points from multidim gaussian." << endl;
	//int n = 0;
	sigma.print();
	Matrix L = sigma.cholesky();	
	//cout << "0,0: " << L.at(0,0) << endl;
	L.print();
	RandomSample rs;
	rs.SetRange(0.,1.);
	vector<double> y;
//	cout << "y size: " << y.size() << endl;
	//for(int n = 0; n < Nsample; n++){
		//y_i ~ N(0,1)
//		cout << "1" << endl;
		y = rs.SampleGaussian(0.,1.,d);
//		cout << "2" << endl;
//		Matrix mat_y(y);
		//cout << "3" << endl;
/*
		//L*y	
		Matrix Ly_mu = Matrix(d,1);
		cout << "4" << endl;
		Ly_mu.mult(L,mat_y);	
		cout << "5" << endl;
		//L*y + mu
		Ly_mu.add(mean);
		cout << "6" << endl;
		//add point to mat
		m_entries.push_back({});
		cout << "7" << endl;
		for(int j = 0; j < d; j++){
			m_entries[n].push_back(Ly_mu.at(n,j));
		} 
		cout << "8" << endl;
//	}

*/
	cout << "SampleNDimGaussian - end" << endl;
}

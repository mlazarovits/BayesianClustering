#include "Matrix.hh"
#include "RandomSample.hh"

Matrix::Matrix(){

}

Matrix::Matrix(int dim){
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
	m_Xdim = x_dim;
	m_Ydim = y_dim;
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

//constructor from inputs

Matrix::~Matrix(){
	for(int i = 0; i < m_dim; i++){
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
}

void Matrix::SetEntry(double val, int i, int j){
	m_entries[i][j] = val;

}

double Matrix::det(){
	if(m_Xdim != m_Ydim){
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


vector<int> Matrix::GetDims(){
	return {m_Xdim, m_Ydim};
}

double Matrix::GetEntry(int i, int j){
	return m_entries[i][j];
}


Matrix Matrix::mult(const Matrix& mat){
	Matrix ret(m_Xdim, dims[1]);
	double val;
	vector<int> dims = mat.GetDims();
	if(dims[0] != m_Ydim){
		cout << "Error: matrices have incompatible dimensions. " << dims[0] << " " m_Ydim << endl;
		return ret;
	}
	for(int i = 0; i < m_Xdim; i++){
		for(int j = 0; j < m_Ydim; j++){
			val = 0;
			for(int k = 0; k < dims[0]; k++){
				val += m_entries[i][k] * mat[k][j];
			}
		}
		ret.SetEntry(val,i,j);
	}
	return ret;
}


Matrix Matrix::transpose(){
	Matrix ret = Matrix(m_Ydim, m_Xdim);
	for(int i = 0; i < m_Xdim; i++){
		for(j = 0; j < m_Ydim; j++){
			ret.SetEntry(m_entries[j][i],i,j);
		}
	}
	return ret;
}
			



}
}


//update entries
//add entries
//clear


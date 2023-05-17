#include "Matrix.hh"
#include "RandomSample.hh"

Matrix::Matrix(){

}

Matrix::Matrix(int dim){
	m_Xdim = dim;
	m_Ydim = dim;
}

Matrix::Matrix(int x_dim, int y_dim){
	m_Xdim = x_dim;
	m_Ydim = y_dim;
}

//constructor from inputs

Matrix::~Matrix(){
	for(int i = 0; i < m_dim; i++){
		m_entries[i].clear();
	}
}

void Matrix::InitRandom(){
	RandomSample rs;
	//TODO: set range by data
	rs.SetRange(0., 100.);
	for(int i = 0; i < m_Xdim; i++){
		vector<double> vec;
		m_entries.push_back(vec);
		for(int j = 0; j < m_Ydim; j++){
			m_entries[i].push_back(rs.SampleFlat());
		}
	}		
}

void Matrix::SetDims(int x_dim, int y_dim){
	m_Xdim = x_dim;
	m_Ydim = y_dim;
}

//update entries
//add entries
//clear


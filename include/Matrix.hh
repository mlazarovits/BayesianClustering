#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
using std::vector;

//inherit from vector of doubles?
class Matrix{
	public:
		//set dimensionality
		Matrix();
		Matrix(int dim);
		Matrix(int x_dim, int y_dim);
		//constructor from inputs
		virtual ~Matrix();

		void InitRandom();	
		void SetDims(int x_dim, int y_dim);
		//update entries
		//add entries
		//clear
	private:
		int m_Xdim;
		int m_Ydim;
		vector<vector<double> m_entries;

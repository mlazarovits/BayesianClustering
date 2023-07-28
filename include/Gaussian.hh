#ifndef GAUSSIAN_HH
#define GAUSSIAN_HH
#include "BasePDF.hh"
#include "Matrix.hh"

class Gaussian : public BasePDF{
	public:
		Gaussian();
		Gaussian(int d);
		Gaussian(Point mu, Matrix cov);
		Gaussian(Matrix mu, Matrix cov);
		virtual ~Gaussian(){ };


		void InitParameters();
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		void UpdateParameters(){ m_mu = m_params["mean"]; m_cov = m_params["cov"]; }	
		double ConjugateEvidence(const Point& x);
		double Prob(const Point& x);
		double Prob(const PointCollection& x);

		Gaussian* mult(BasePDF* p1){ 
			//check that p1 is a Gaussian by looking at parameters
			Matrix mu2 = p1->GetParameter("mean");
			Matrix mu1 = m_mu;

			Matrix cov2 = p1->GetParameter("cov");
			cov2.invert(cov2);
			Matrix cov1 = Matrix(m_dim, m_dim);
			cov1.invert(m_cov);
			

			Gaussian* ret = new Gaussian(m_dim);
			if(mu2.empty() || cov2.empty()){ cout << "Error: input must be same as derived class (Gaussian)" << endl; return ret; }

			//form from: https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf (Ch. 8.1.8)
			Matrix new_cov = Matrix(m_dim, m_dim);
			new_cov.add(cov1,cov2);
			
			Matrix new_mu1 = Matrix(m_dim, 1);
			Matrix new_mu2 = Matrix(m_dim, 1);
			Matrix new_mu = Matrix(m_dim, 1);
			
			new_mu1.mult(cov1,mu1);
			new_mu2.mult(cov2,mu2);

			new_mu.add(new_mu1, new_mu2);
			new_mu.mult(new_cov,new_mu);

			//use added covariances in this distribution
			Gaussian* gaus_coeff = new Gaussian(mu2, new_cov);
			//then invert added covariances for overall distribution
			new_cov.invert(new_cov); 

			Point x = mu1.MatToPoints().at(0);
			double coeff = gaus_coeff->Prob(x);	

			ret->SetParameter("mean",new_mu);
			ret->SetParameter("cov",new_cov);
			ret->SetPrefactor(coeff);
	
			return ret;	
		};


	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif

#ifndef VARGAUSSIANMIXTURE_HH
#define VARGAUSSIANMIXTURE_HH

#include "BasePDFMixture.hh"

class VarGaussianMixture : public BasePDFMixture{

	public:
		VarGaussianMixture();
		virtual ~VarGaussianMixture();

		void Initialize(unsigned long long seed = 111);

		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();









};


#endif






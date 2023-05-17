#ifndef BASEPDFMIXTURE_HH
#define BASEPDFMIXTURE_HH


class BasePDFMixture{
	public:
		BasePDFMixture();
		
		virtual ~BasePDFMixture(){ };
		
		virtual void Initialize() = 0;
		//E-step
		virtual void CalculatePosterior() = 0;
		//M-step
		virtual void UpdateParameters(BasePDFMixture* post) = 0;
		//some stuff for BHC













#endif

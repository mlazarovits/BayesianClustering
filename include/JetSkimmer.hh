#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "TFile.h"
#include "BaseSkimmer.hh"
#include "BaseProducer.hh"


class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer();
		virtual ~JetSkimmer();

		//get rechits from file to cluster
		JetSkimmer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		
		void CleaningSkim(){ };
		void Skim();

		private:
			BaseProducer* _prod;






};
#endif

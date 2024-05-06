#ifndef SampleWeight_HH
#define SampleWeight_HH

#include <string>
#include <map>
#include <TFile.h>

using std::string;
using std::map;

class SampleWeight{
	public:
		SampleWeight();
		void Init();
		virtual ~SampleWeight(){ };
		struct weights{
			double scale;
			double xsec;
		
			weights(){
				scale = 0.;
				xsec = 0.;
			}
			void SetScale(double s){scale = s;}
			void SetXsec(double s){xsec = s; }
		};

		double GetScale(string samp){ return _sampleToWeights[samp].scale; }
		double GetXsec(string samp){ return _sampleToWeights[samp].xsec; }

		double GetScale(TFile* f){
			string name = f->GetName();
			string key;
			for(auto it = _sampleToWeights.begin(); it != _sampleToWeights.end(); it++){
				key = it->first;
				if(name.find(key) != string::npos) return it->second.scale;
			}
			return 1.;
		}
		
		double GetXsec(TFile* f){
			string name = f->GetName();
			string key;
			for(auto it = _sampleToWeights.begin(); it != _sampleToWeights.end(); it++){
				key = it->first;
				if(name.find(key) != string::npos) return it->second.xsec;
			}
			return 1.;
		}
		
		void GetWeights(TFile* f, double& scale, double& xsec){
			string name = f->GetName();
			string key;
			for(auto it = _sampleToWeights.begin(); it != _sampleToWeights.end(); it++){
				key = it->first;
				if(name.find(key) != string::npos){
					scale = it->second.scale;
					xsec = it->second.xsec;
					return;
				} 
			}
			scale = 1.;
			xsec = 1.;
		}

		weights GetWeights(string samp){return _sampleToWeights[samp]; }

			

	private:
		map<string,weights> _sampleToWeights;

};
#endif

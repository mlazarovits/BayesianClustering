#include "SampleWeight.hh"

SampleWeight::SampleWeight(){
	//make map
	_sampleToWeights["QCD_HT50to100"] = weights();
	_sampleToWeights["QCD_HT100to200"] = weights();
	_sampleToWeights["QCD_HT200to300"] = weights();
	_sampleToWeights["QCD_HT300to500"] = weights();
	_sampleToWeights["QCD_HT500to700"] = weights();
	_sampleToWeights["QCD_HT700to1000"] = weights();
	_sampleToWeights["QCD_HT1000to1500"] = weights();
	_sampleToWeights["QCD_HT1500to2000"] = weights();
	_sampleToWeights["QCD_HT2000toInf"] = weights();
	
}


void SampleWeight::Init(){
	//scale
	_sampleToWeights["QCD_HT50to100"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT100to200"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT200to300"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT300to500"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT500to700"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT700to1000"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT1000to1500"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT1500to2000"].SetLumiScale(1.);
	_sampleToWeights["QCD_HT2000toInf"].SetLumiScale(1.);


	//xsec
	_sampleToWeights["QCD_HT50to100"].SetXsec(187700000.0);
	_sampleToWeights["QCD_HT100to200"].SetXsec(23670000.0);
	_sampleToWeights["QCD_HT200to300"].SetXsec(1554000.0);
	_sampleToWeights["QCD_HT300to500"].SetXsec(323800.0);
	_sampleToWeights["QCD_HT500to700"].SetXsec(30010.0);
	_sampleToWeights["QCD_HT700to1000"].SetXsec(6352.0);
	_sampleToWeights["QCD_HT1000to1500"].SetXsec(1096.0);
	_sampleToWeights["QCD_HT1500to2000"].SetXsec(99.12);
	_sampleToWeights["QCD_HT2000toInf"].SetXsec(20.20);



}


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

	_sampleToWeights["GMSB_L-100TeV"] = weights();  
	_sampleToWeights["GMSB_L-150TeV"] = weights();  
	_sampleToWeights["GMSB_L-200TeV"] = weights();  
	_sampleToWeights["GMSB_L-250TeV"] = weights();  
	_sampleToWeights["GMSB_L-300TeV"] = weights();  
	_sampleToWeights["GMSB_L-350TeV"] = weights();  
	_sampleToWeights["GMSB_L-400TeV"] = weights();  
	_sampleToWeights["GMSB_L-500TeV"] = weights();  
	_sampleToWeights["GMSB_L-600TeV"] = weights();  
	
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

	_sampleToWeights["GMSB_L-100TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-150TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-200TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-250TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-300TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-350TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-400TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-500TeV"].SetLumiScale(1.);  
	_sampleToWeights["GMSB_L-600TeV"].SetLumiScale(1.);  
	
	_sampleToWeights["GJets_HT40to100"].SetLumiScale(1.);
	_sampleToWeights["GJets_HT100to200"].SetLumiScale(1.);
	_sampleToWeights["GJets_HT200to400"].SetLumiScale(1.);
	_sampleToWeights["GJets_HT400to600"].SetLumiScale(1.);
	_sampleToWeights["GJets_HT600toInf"].SetLumiScale(1.);

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

	_sampleToWeights["GMSB_L-100TeV"].SetXsec(2.16);  
	_sampleToWeights["GMSB_L-150TeV"].SetXsec(0.228);  
	_sampleToWeights["GMSB_L-200TeV"].SetXsec(0.0445);  
	_sampleToWeights["GMSB_L-250TeV"].SetXsec(0.0126);  
	_sampleToWeights["GMSB_L-300TeV"].SetXsec(0.00445);  
	_sampleToWeights["GMSB_L-350TeV"].SetXsec(0.00178);  
	_sampleToWeights["GMSB_L-400TeV"].SetXsec(0.000778);  
	_sampleToWeights["GMSB_L-500TeV"].SetXsec(0.000165);  
	_sampleToWeights["GMSB_L-600TeV"].SetXsec(0.0);  

	_sampleToWeights["GJets_HT40to100"].SetXsec(18620.0);
	_sampleToWeights["GJets_HT100to200"].SetXsec(8625.0);
	_sampleToWeights["GJets_HT200to400"].SetXsec(2196.0);
	_sampleToWeights["GJets_HT400to600"].SetXsec(258.0);
	_sampleToWeights["GJets_HT600toInf"].SetXsec(85.18);


}

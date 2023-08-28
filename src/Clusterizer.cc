#include "Clusterizer.hh"
#include "VarClusterViz3D.hh"
#include "BayesHierCluster.hh"
#include "FullViz3D.hh"
#include "KMeansCluster.hh"
#include <TSystem.h>

Clusterizer::Clusterizer(){
	m_nJets = 0;
};

Clusterizer::Clusterizer(vector<Jet> jets){
	//add jet info to pointcollection
	//match jets to points by index
	m_oldJets = jets;
}


Clusterizer::~Clusterizer(){ }

void Clusterizer::Cluster(Jet jet, double alpha, double thresh, bool viz, int verb, string fname){
	PointCollection* points = new PointCollection();
	jet.GetEtaPhiConstituents(*points);
	//Bayesian Hierarchical Clustering algo
	BayesHierCluster* bhc = new BayesHierCluster(alpha);
	bhc->SetThresh(thresh);
	bhc->AddData(points);
	bhc->SetVerbosity(verb);
	if(!_params.empty()) bhc->SetPriorParameters(_params);
	if(!_data_smear.empty()) bhc->SetDataSmear(_data_smear);
	//each node is a jet - a mixture of gaussians (subjets)
	vector<node*> tree = bhc->Cluster();
	if(viz){
		if(fname.empty()) fname = "plots/test";
		FullViz3D plots = FullViz3D(tree);
		plots.SetVerbosity(verb);
		plots.Write(fname);
	}
	//get parameters from subjets - GMM components
	if(verb > 0) cout << tree.size() << " jets found." << endl;
	for(int n = 0; n < (int)tree.size(); n++){
		BasePDFMixture* model = tree[n]->model;
		if(verb > 1){
			cout << model->GetNClusters() << " subjets found in jet " << n << endl;
			cout << "Estimated prior parameters" << endl;
			map<string, Matrix> params;
			for(int i = 0; i < model->GetNClusters(); i++){
				params = model->GetPriorParameters(i);	
				cout << "weight " << i << ": " << params["pi"].at(0,0) << " alpha " << params["alpha"].at(0,0) << endl;
				cout << "mean " << i << endl;
				params["m"].Print();
				cout << "Gaus scale " << params["scale"].at(0,0) << endl;
				cout << "dof " << params["dof"].at(0,0) << endl;
				cout << "scalemat " << i << endl;
				params["scalemat"].Print();
				cout << "mean " << i << endl;
				params["mean"].Print();
				cout << "cov " << i << endl;
				params["cov"].Print();
				params.clear();
			}
		}
	}

}




//crack open Jet and get underlying points
void Clusterizer::FindSubjets(Jet jet, double alpha, double thresh, bool viz, int verb, int maxK, string fname){
	//create GMM model
	PointCollection* points = new PointCollection();
	jet.GetEtaPhiConstituents(*points);

	vector<double> weights;
	jet.GetEnergies(weights);
	
	points->SetWeights(weights);
	fname += "_Eweighted";

	GaussianMixture* gmm = new GaussianMixture(maxK);
	
	gmm->SetData(points);
	gmm->SetAlpha(alpha);
	gmm->InitParameters();
	gmm->InitPriorParameters();


	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,maxK);
	algo->SetThresh(thresh);
	

	map<string, vector<Matrix>> params;
	
	if(fname.empty()) fname = "plots/test";
	if(gSystem->AccessPathName((fname).c_str())){
		gSystem->Exec(("mkdir -p "+fname).c_str());
	}
	else{
		gSystem->Exec(("rm -rf "+fname).c_str());
		gSystem->Exec(("mkdir -p "+fname).c_str());

	}
	VarClusterViz3D cv3D;
	if(viz){ cv3D = VarClusterViz3D(algo);
		cv3D.UpdatePosterior();
		cv3D.WriteJson(fname+"/it0");
		}
	//loop
	double dLogL, newLogL;
	double LogLthresh = 0.01;
	double oldLogL = algo->EvalLogL();
	////////run EM algo////////
	//maximum of 50 iterations
	for(int it = 0; it < 50; it++){
	
		//E step
		algo->Estimate();
		//M step
		algo->Update();
		
		//Plot
		if(viz){
			cv3D.UpdatePosterior();
			cv3D.WriteJson(fname+"/it"+std::to_string(it+1));
		}
		//Check for convergence
		newLogL = algo->EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << endl;
			return;
		}
		dLogL = oldLogL - newLogL;
		if(viz) cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(fabs(dLogL) < LogLthresh){
			if(viz){
				cout << "Reached convergence at iteration " << it+1 << endl;
			}
			break;
		}
		oldLogL = newLogL;
	}
	if(viz){
		cout << "Estimated parameters" << endl;
		map<string, Matrix> params;
		for(int i = 0; i < gmm->GetNClusters(); i++){
			params = gmm->GetPriorParameters(i);	
			cout << "weight " << i << ": " << params["pi"].at(0,0) << endl;
			cout << "mean " << i << endl;
			params["mean"].Print();
			cout << "cov " << i << endl;
			params["cov"].Print();
			params.clear();
		}

	}
}











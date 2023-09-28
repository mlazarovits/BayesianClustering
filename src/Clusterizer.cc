#include "Clusterizer.hh"
#include "VarClusterViz3D.hh"
#include "BayesHierCluster.hh"
#include "FullViz3D.hh"
#include "KMeansCluster.hh"
#include <TSystem.h>

Clusterizer::Clusterizer(){
	m_nJets = 0;
	_bhcAlpha = 0;
	_emAlpha = 0;
	_thresh = 0;
	_verb = 0;
	_maxK = 0;
	_weighted = false;
	_smeared = false;
	_distconst = false;
};

Clusterizer::Clusterizer(vector<Jet> jets){
	//add jet info to pointcollection
	//match jets to points by index
	m_oldJets = jets;
	_bhcAlpha = 0;
	_emAlpha = 0;
	_thresh = 0;
	_verb = 0;
	_maxK = 0;
	_weighted = false;
	_smeared = false;
	_distconst = false;
}


Clusterizer::~Clusterizer(){ }

vector<node*> Clusterizer::Cluster(Jet jet, string fname){
	_points = new PointCollection();
	jet.GetEtaPhiConstituents(*_points);

	
	if(_weighted){
		vector<double> weights;
		jet.GetEnergies(weights);
		//need to transfer from GeV (energy) -> unitless (number of points)
		double gev = 0;
		for(int i = 0; i < (int)weights.size(); i++) gev += weights[i];
		gev = gev/(double)weights.size(); //k = sum_n E_n/n pts
		for(int i = 0; i < (int)weights.size(); i++) weights[i] /= gev; //sums to n pts, w_n = E_n/k  
		_points->SetWeights(weights);

	}

	//Bayesian Hierarchical Clustering algo
	BayesHierCluster* bhc = new BayesHierCluster(_bhcAlpha);
	bhc->SetSubclusterAlpha(_emAlpha);	
	bhc->SetVerbosity(_verb);
	//set configs
	if(_smeared) bhc->SetDataSmear(_data_smear);
	bhc->SetThresh(_thresh);
	bhc->AddData(_points);
	if(!_params.empty()) bhc->SetPriorParameters(_params);
	
	if(_distconst){
		//int d = dimension, double c = threshold, double a = lower bound, double b = upper bound
		//setting constraint of pi/2 on phi -> dphi must be at least pi/2 between two clusters
		bhc->SetDistanceConstraint(0., acos(-1)/2);
	}
	//run algo
	//each node is a jet - a mixture of gaussians (subjets)
	vector<node*> tree = bhc->Cluster();


	//plotting
	if(!fname.empty()){
		FullViz3D plots = FullViz3D(tree);
		plots.SetVerbosity(_verb);
		plots.Write(fname);
		//plot individual clusters in local coordinates here
		//for(int k = 0; k < tree; k++)
		//plots.SetLocalCoords(true);
		//json object root = plots.WriteNode(tree[i]);
		//std::ofstream file;
		//file.open(filename+to_string(k)+".json");
		////sets 4 space indent
		//file << std::setw(4) << root << endl;
		//cout << "Writing to: " << filename << ".json" << endl;
	}


	//get parameters from subjets - GMM components
	if(_verb > 0) cout << tree.size() << " jets found." << endl;
	for(int n = 0; n < (int)tree.size(); n++){
		BasePDFMixture* model = tree[n]->model;
		if(_verb > 3){
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

	return tree;

}




//crack open Jet and get underlying points
GaussianMixture* Clusterizer::FindSubjets(Jet jet, string fname){
	//create GMM model
	_points = new PointCollection();
	jet.GetEtaPhiConstituents(*_points);

	if(_weighted){
		vector<double> weights;
		jet.GetEnergies(weights);
		//need to transfer from GeV (energy) -> unitless (number of points)
		double gev = 0;
		for(int i = 0; i < (int)weights.size(); i++) gev += weights[i];
		gev = gev/(double)weights.size(); //k = sum_n E_n/n pts
		for(int i = 0; i < (int)weights.size(); i++) weights[i] /= gev; //sums to n pts 
		_points->SetWeights(weights);
	}

	GaussianMixture* gmm = new GaussianMixture(_maxK);
	
	gmm->SetData(_points);
	gmm->SetAlpha(_emAlpha);
	gmm->SetVerbosity(_verb);
	gmm->InitParameters();
	gmm->InitPriorParameters();

	if(_smeared){
		gmm->SetDataSmear(_data_smear);
	}

	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,_maxK);
	algo->SetThresh(_thresh);
	

	map<string, vector<Matrix>> params;
	bool viz = false;
	if(!fname.empty()){
		viz = true;
	}

	VarClusterViz3D cv3D;
	if(viz){ cv3D = VarClusterViz3D(algo);
		cv3D.SetVerbosity(_verb);
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
			return gmm;
		}
		dLogL = oldLogL - newLogL;
		if(_verb > 0) cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(fabs(dLogL) < LogLthresh){// || dLogL > 0){
			if(_verb > 0){
				cout << "Reached convergence at iteration " << it+1 << endl;
			}
			break;
		}
		oldLogL = newLogL;
	}
	if(_verb > 1){
		cout << "Estimated parameters" << endl;
		map<string, Matrix> params;
		for(int i = 0; i < gmm->GetNClusters(); i++){
			params = gmm->GetPriorParameters(i);	
			cout << "weight " << i << ": " << params["pi"].at(0,0) << " alpha: " << params["alpha"].at(0,0) << endl;
			cout << "mean " << i << endl;
			params["mean"].Print();
			cout << "cov " << i << endl;
			params["cov"].Print();
			params.clear();
		}

	}

	return gmm;
}











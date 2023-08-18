#include "JetClusterizer.hh"
#include "VarClusterViz3D.hh"
#include "BayesHierCluster.hh"
#include "FullViz3D.hh"
#include <TSystem.h>

JetClusterizer::JetClusterizer(){
	m_nJets = 0;
};

JetClusterizer::JetClusterizer(vector<Jet> jets){
	//add jet info to pointcollection
	//match jets to points by index
	m_oldJets = jets;
}


JetClusterizer::~JetClusterizer(){ }

void JetClusterizer::Cluster(Jet jet, double alpha, double thresh, bool viz, int verb){
	PointCollection* points = new PointCollection();
	jet.GetEtaPhiConstituents(*points);
	//Bayesian Hierarchical Clustering algo
	BayesHierCluster* bhc = new BayesHierCluster(alpha);
	bhc->SetThresh(thresh);
	bhc->AddData(points);
	bhc->SetVerbosity(verb);
	//each node is a jet - a mixture of gaussians (subjets)
	vector<node*> tree = bhc->Cluster();
	if(viz){
		string fname = "plots/jettest";
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
			cout << "Estimated subjet parameters" << endl;
			map<string, Matrix> params;
			for(int i = 0; i < model->GetNClusters(); i++){
				params = model->GetPriorParameters(i);	
				cout << "weight " << i << ": " << params["pi"].at(0,0) << endl;
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
GaussianMixture* JetClusterizer::FindSubjets(PointCollection* points, double thresh, int maxNit, int maxK, bool viz, double a, PointCollection* seeds){

	//create GMM model
	GaussianMixture* gmm = new GaussianMixture(maxK);
	
	Point norm_scale = points->Normalize();

	gmm->SetData(points);
	gmm->SetAlpha(a);
	if(seeds != nullptr) gmm->InitParameters(*seeds);
	else gmm->InitParameters();
	gmm->InitPriorParameters();


	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,maxK);
	algo->SetThresh(thresh);
	

	map<string, vector<Matrix>> params;
	
	string fname = "plots/jetTest/"; 
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
		cv3D.WriteJson(fname+"it0");
		}
	//loop
	double dLogL, newLogL;
	double LogLthresh = 0.01;
	double oldLogL = algo->EvalLogL();
	////////run EM algo////////
	for(int it = 0; it < maxNit; it++){
	
		//E step
		algo->Estimate();
		//M step
		algo->Update();
		
		//Plot
		if(viz){
			cv3D.UpdatePosterior();
			cv3D.WriteJson(fname+"it"+std::to_string(it+1));
		}
		//Check for convergence
		newLogL = algo->EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << endl;
			return gmm;
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

	return gmm;
}



//crack open Jet and get underlying points
vector<Jet> JetClusterizer::FindSubjets_etaPhi(Jet jet, double thresh, int maxNit, int maxK, bool viz, double a){
	vector<Jet> subjets;
	//vector<JetPoint> rhs;
	//jet.GetConstituents(rhs);
	//Point vtx = jet.GetVertex();

	//initialize vector of subjets
	PointCollection* points = new PointCollection();
	jet.GetEtaPhiConstituents(*points);
	int n_pts = points->GetNPoints();

	/*	
	//find maxK points with biggest energy - 4th dim
	points.Sort(3);
	PointCollection seeds;
	PointCollection newpts;
	vector<double> pt;
	for(int k = 0; k < points.GetNPoints(); k++){
		//discard energy dimension
		pt.push_back(points.at(k).at(0));
		pt.push_back(points.at(k).at(1));
		pt.push_back(points.at(k).at(2));

		if(k >= points.GetNPoints() - maxK) seeds += Point(pt);
		newpts += Point(pt);
		pt.clear();
	}
	points.Clear();
	*/
	GaussianMixture* gmm = FindSubjets(points, thresh, maxNit, maxK, viz, a);
	int nsubjets = gmm->GetNClusters();


	return subjets;

}

//crack open Jet and get underlying points
vector<Jet> JetClusterizer::FindSubjets_XYZ(Jet jet, double thresh, int maxNit, int maxK, bool viz, double a){
	vector<Jet> subjets;

	//initialize vector of subjets
	PointCollection* points = new PointCollection();
	jet.GetXYZConstituents(*points);
	int n_pts = points->GetNPoints();

	/*	
	//find maxK points with biggest energy - 4th dim
	points.Sort(3);
	PointCollection seeds;
	PointCollection newpts;
	vector<double> pt;
	for(int k = 0; k < points.GetNPoints(); k++){
		//discard energy dimension
		pt.push_back(points.at(k).at(0));
		pt.push_back(points.at(k).at(1));
		pt.push_back(points.at(k).at(2));

		if(k >= points.GetNPoints() - maxK) seeds += Point(pt);
		newpts += Point(pt);
		pt.clear();
	}
	points.Clear();
	*/
	GaussianMixture* gmm = FindSubjets(points, thresh, maxNit, maxK, viz, a);
	int nsubjets = gmm->GetNClusters();


	return subjets;
}







#include "JetClusterizer.hh"
#include "VarClusterViz3D.hh"

JetClusterizer::JetClusterizer(){
	m_nJets = 0;
};

JetClusterizer::JetClusterizer(vector<Jet> jets){
	//add jet info to pointcollection
	//match jets to points by index
	m_oldJets = jets;
}


JetClusterizer::~JetClusterizer(){ }

void JetClusterizer::Cluster(){
	/*
	if(m_nJets > 1){
		//BHC algo
		m_newJets = MergeJets(m_oldJets);
		vector<Jet> subJets;
		int p1, p2;
		for(int j = 0; j < (int)m_newJets.size(); j++){
			subJets = FindSubjets(m_newJets[j]);
			//set parents for merge (newJet) j
			BHC.GetSubTrees(p1, p2, j);
			m_newJets[j].SetParents(m_oldJets[p1],m_oldJets[p2]);
		}	
		m_nJets = (int)newJets.size();
		//TODO: make sure oldJets have children set as newJets in tree
		m_tree.AddLayer(m_newJets);
		m_oldJets.clear();
		m_oldJets = m_newJets;
		m_newJets.clear();
		subJets.clear();
		Cluster();
	}
	else{
		return;
	}


*/

}

//crack open Jet and get underlying points
vector<Jet> JetClusterizer::FindSubjets_etaPhi(Jet jet, double LogLthresh, int maxNit, int maxK, bool viz){
	//initialize vector of subjets
	vector<Jet> subjets;
	vector<JetPoint> rhs;
	Point vtx = jet.GetVertex();

	PointCollection points;
	VarGaussianMixture vgmm(maxK);
	jet.GetEtaPhiConstituents(points);
	jet.GetConstituents(rhs);
	int n_pts = points.GetNPoints();
cout << n_pts << " constituents to cluster" << endl;	
	
	vgmm.AddData(points);
	vgmm.Initialize();
	VarClusterViz3D cv3D;
	string fname = "jetTest"; 
	if(viz) cv3D = VarClusterViz3D(&vgmm, fname);
	//loop
	double dLogL, newLogL, oldLogL;
	cout << "maxNit: " << maxNit << endl;
	////////run EM algo////////
	for(int it = 0; it < maxNit; it++){
		oldLogL = vgmm.EvalLogL();
		
		//E step
		vgmm.CalculatePosterior();
		//M step
		vgmm.UpdateParameters();
		
		//Plot
		if(viz){
			cv3D.UpdatePosterior();
			cv3D.AddAnimationPlots("it"+std::to_string(it));
		}
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return subjets;
		}
		dLogL = fabs(oldLogL - newLogL);
		if(viz) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(dLogL < LogLthresh){
			if(viz)
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
	}
	if(viz)	cv3D.Write();
	int nsubjets = vgmm.GetNClusters(LogLthresh*10.);

	//TODO: consider edge case where r_nk = r_nk' for all k == k' (all k entries for a point are equal)
	//assign points to found subjets (clusters)
	for(int k = 0; k < nsubjets; k++){
		subjets.push_back(Jet());
		subjets[k].SetVertex(vtx);
	}
	int k_assign;
	for(int n = 0; n < n_pts; n++){
		//find k where post_nk is max for point n:
		k_assign = vgmm.GetMaxPointAssignment(n);
		//assign point n to subjet k - subjets[postMax_k] += rh[n];
		subjets[k_assign].add(rhs[n]);
	}
	
	return subjets;

}

//crack open Jet and get underlying points
vector<Jet> JetClusterizer::FindSubjets_XYZ(Jet jet, double LogLthresh, int maxNit, int maxK, bool viz){
	//initialize vector of subjets
	vector<Jet> subjets;
	vector<JetPoint> rhs;

	PointCollection points;
	VarGaussianMixture vgmm(maxK);
	jet.GetPointConstituents(points);
	jet.GetConstituents(rhs);
	int n_pts = points.GetNPoints();
cout << n_pts << " constituents to cluster" << endl;	
	
	vgmm.AddData(points);
	vgmm.Initialize();
	VarClusterViz3D cv3D;
	string fname = "jetTest"; 
	if(viz) cv3D = VarClusterViz3D(&vgmm, fname);
	//loop
	double dLogL, newLogL, oldLogL;
	////////run EM algo////////
	for(int it = 0; it < maxNit; it++){
		oldLogL = vgmm.EvalLogL();
		
		//E step
		vgmm.CalculatePosterior();
		//M step
		vgmm.UpdateParameters();
		
		//Plot
		if(viz){
			cv3D.UpdatePosterior();
			cv3D.AddAnimationPlots("it"+std::to_string(it));
		}
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return subjets;
		}
		dLogL = fabs(oldLogL - newLogL);
		if(viz) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(dLogL < LogLthresh){
			if(viz)
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
	}
	if(viz)	cv3D.Write();
	int nsubjets = vgmm.GetNClusters(LogLthresh*10.);

	//TODO: consider edge case where r_nk = r_nk' for all k == k' (all k entries for a point are equal)
	//assign points to found subjets (clusters)
	for(int k = 0; k < nsubjets; k++){
		subjets.push_back(Jet());
	}
	int k_assign;
	for(int n = 0; n < n_pts; n++){
		//find k where post_nk is max for point n:
		k_assign = vgmm.GetMaxPointAssignment(n);
		//assign point n to subjet k - subjets[postMax_k] += rh[n];
		subjets[k_assign].add(rhs[n]);
	}
	
	return subjets;

}







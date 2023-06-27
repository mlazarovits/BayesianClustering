#include "JetClusterizer.hh"
#include "VarClusterViz2D.hh"

JetClusterizer::JetClusterizer(){
	m_maxK = 0;
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
vector<Jet> JetClusterizer::FindSubjets(Jet jet, double LogLthresh, int maxNit, bool viz){
	//initialize vector of subjets
	vector<Jet> subjets;

	PointCollection points;
	VarGaussianMixture vgmm(m_maxK);
	//match points to jets by idx
	jet.GetConstituents(points);
	int n_pts = points.GetNPoints();
	
	vgmm.AddData(points);
	vgmm.Initialize();
	VarClusterViz2D cv2D;
	string fname = "test"; 
	if(viz) cv2D = VarClusterViz2D(&vgmm, fname);
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
			cv2D.UpdatePosterior();
			cv2D.AddPlot("it"+std::to_string(it));
		}
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return subjets;
		}
		dLogL = fabs(oldLogL - newLogL);
	//	cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
	//	if(dLogL < LogLThresh){
	//		cout << "Reached convergence at iteration " << it << endl;
	//		break;
	//	}
		
	}
	if(viz)	cv2D.Write();
	int nsubjets = vgmm.GetNClusters(LogLthresh*10.);

	//TODO: consider edge case where r_nk = r_nk' for all k == k' (all k entries for a point are equal)
	//assign points to found subjets (clusters)
	for(int k = 0; k < nsubjets; k++)
		subjets.push_back(Jet());
	//for(int n = 0; n < m_n; 
	

	return subjets;
}







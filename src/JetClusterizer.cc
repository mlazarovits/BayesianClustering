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

	//normalize points
	//find max/min in each dim
	double etaMin = points.min(0);
	double etaMax = points.max(0);
	double phiMin = points.min(1);
	double phiMax = points.max(1);
	double tMin = points.min(2);
	double tMax = points.max(2);
	
cout << "tMin: " << tMin << " tMax: " << tMax << endl;

	Matrix matPts = Matrix(3, points.GetNPoints());
	matPts.PointsToMat(points);
	//translate by -xmin - lower bound to 0
	Matrix minTrans = Matrix(3,points.GetNPoints());	
	for(int i = 0; i < points.GetNPoints(); i++){
		minTrans.SetEntry(-etaMin,0,i);
		minTrans.SetEntry(-phiMin,1,i);
		minTrans.SetEntry(-tMin,2,i);
	}

	//scale by xmax-xmin - normalize to 1
	Matrix normScale = Matrix(3,3);	
	normScale.InitEmpty();
	normScale.SetEntry(1./(etaMax-etaMin),0,0);
	normScale.SetEntry(1./(phiMax-phiMin),1,1);
	normScale.SetEntry(1./(tMax-tMin),2,2);
	//translate first
	matPts.add(minTrans);
	
	//then scale
	matPts.mult(normScale,matPts);

	PointCollection transfPts = matPts.MatToPoints();

	//cout << "points" << endl;
	//points.Print();

	//cout << "transfPoints" << endl;
	//transfPts.Print();

	//vgmm.AddData(points);
	vgmm.AddData(transfPts);
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
			cv3D.AddAnimation("it"+std::to_string(it));
		}
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return subjets;
		}
		dLogL = fabs(oldLogL - newLogL);
		if(viz) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << "\n" << endl;
		if(dLogL < LogLthresh){
			if(viz){
				cout << "Reached convergence at iteration " << it << endl;
			}
			break;
		}
	}
	//if(viz)	cv3D.Write();
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
			cv3D.AddAnimation("it"+std::to_string(it));
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
	//if(viz)	cv3D.Write();
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







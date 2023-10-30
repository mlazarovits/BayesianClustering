#include "FullViz3D.hh"

FullViz3D::FullViz3D(const vector<node*>& nodes){
	for(int i = 0; i < (int)nodes.size(); i++){
		if(nodes[i] == nullptr) continue;
		_nodes.push_back(nullptr);
		node* x = new node(*(nodes[i]));
		_nodes[(int)_nodes.size() - 1] = nodes[i];
	}
		
}



json FullViz3D::WriteNode(node* n){
	json subcluster = json::object();
	json subclusters = json::object();
	json cluster = json::object();
	json data = json::object();

	vector<double> color, x, y, z, w;
	
	vector<double> eigenVec_0, eigenVec_1, eigenVec_2;
	
	BasePDFMixture* pdfmodel = n->model;
	int kmax = pdfmodel->GetNClusters();
	PointCollection* points = n->points;
	if(points->GetNPoints() == 0){
		return cluster;
	}

	for(int i = 0; i < points->GetNPoints(); i++){
		//eta
		x.push_back(points->at(i).Value(0));
		//phi
		y.push_back(points->at(i).Value(1));
		//time
		z.push_back(points->at(i).Value(2));
		//weight - untransfererd (in GeV (sum_n E_n*r_nk))
		w.push_back(points->at(i).Weight()*_transf);

	}
	data["x"] = x;
	data["y"] = y;
	data["z"] = z;
	data["w"] = w;
	cluster["data"] = data;

	//set coords for parameter circles
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	vector<double> cnts;
	pdfmodel->GetNorms(cnts);
	double pisum = 0;
	
	double x0, y0, z0;	
	map<string, Matrix> cluster_params;
	for(int k = 0; k < kmax; k++){
		cluster_params = pdfmodel->GetPriorParameters(k);
		x0 = cluster_params["mean"].at(0,0);
		y0 = cluster_params["mean"].at(1,0);
		z0 = cluster_params["mean"].at(2,0);

		cluster_params["cov"].eigenCalc(eigenVals, eigenVecs);
		for(int i = 0; i < 3; i++){
			eigenVec_0.push_back(eigenVecs[0].at(i,0));
			eigenVec_1.push_back(eigenVecs[1].at(i,0));
			eigenVec_2.push_back(eigenVecs[2].at(i,0));
		}

	//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
		pisum += cluster_params["pi"].at(0,0);
		subcluster["mixing_coeff"] = cluster_params["pi"].at(0,0);
		subcluster["mu_x"] = x0;
		subcluster["mu_y"] = y0;
		subcluster["mu_z"] = z0;

		subcluster["eigenVal_0"] = eigenVals[0];	
		subcluster["eigenVal_1"] = eigenVals[1];	
		subcluster["eigenVal_2"] = eigenVals[2];	
		
		subcluster["eigenVec_0"] = eigenVec_0;	
		subcluster["eigenVec_1"] = eigenVec_1;	
		subcluster["eigenVec_2"] = eigenVec_2;	
		
		//color for subcluster will be total energy (sum_n E_n*r_nk)
		subcluster["color"] = _transf*cnts[k];
	
		subclusters["subcluster_"+std::to_string(k)] = subcluster;


		eigenVec_0.clear();
		eigenVec_1.clear();
		eigenVec_2.clear();
	}

	cluster["subclusters"] = subclusters;
	return cluster;

}



void FullViz3D::orderTree(node* n, int level, map<int, NodeStack> &map){
	//either head or end node
	if(n->val == -1) return;
	//a node that's been deleted
	if(n == nullptr) return;
	//a mirror node
	if(n->points->mean().at(1) < 0.0 || n->points->mean().at(1) >= 2*acos(-1)) return;
	map[level].push(n);


	orderTree(n->l, level + 1, map);
	orderTree(n->r, level + 1, map);

}

//level{ tree_0{ }, tree_1{ }, ...}
json FullViz3D::WriteLevels(){
	json levels = json::object();
	json trees = json::object();
	json clusters = json::object();

	int nTrees = 0;
	
	//create a node - level map for each tree
	vector<map<int, NodeStack>> tree_maps;
	for(int i = 0; i < (int)_nodes.size(); i++){
		map<int, NodeStack> tree_map;
		if(_nodes[i] == nullptr){ continue; }
		//a mirror node
		if(_nodes[i]->points->mean().at(1) < 0.0 || _nodes[i]->points->mean().at(1) >= 2*acos(-1)) continue;
		orderTree(_nodes[i], 0, tree_map);
		tree_maps.push_back(tree_map);
		nTrees++;	
	}
	if(tree_maps.empty()) return levels;

    auto pr = std::max_element(tree_maps.begin(), tree_maps.end(), [](const auto &x, const auto &y) {
                    return x.rbegin()->first < y.rbegin()->first;
                });
	int nLevels = pr->rbegin()->first;
if(_verb > 1) cout << "max: " << nLevels << " levels with " << nTrees << " trees" << endl;
	//write a json for each level per tree
	for(int l = 0; l < nLevels+1; l++){
		if(_verb > 2) cout << "Level " << l << ": " << endl;
		int npts = 0;
		for(int t = 0; t < nTrees; t++){
			//only write if level exists in tree
			if(l <= tree_maps[t].rbegin()->first){
				if(_verb > 2) cout << "  Tree " << t << ": " << endl;
				int j = 0;
				//loop through nodes (clusters) at this level for this tree
				//can't set to empty because need to evaluate last node popped off below
				while(!tree_maps[t][l].empty()){
					node* n = tree_maps[t][l].pop();
					if(_verb > 2) cout << "    node " << j << " - number of points: " << n->points->GetNPoints() << endl; 
					//a mirror node
					if(n->points->mean().at(1) < 0.0 || n->points->mean().at(1) >= 2*acos(-1)) continue;
					//if there is a cluster with one point in tree_maps[t][l] (a leaf) add it to tree_maps[t][l+1]
					if(n->points->GetNPoints() == 1 && l <= tree_maps[t].rbegin()->first){
						tree_maps[t][l+1].push(n);
					}
					npts += n->points->GetNPoints();
					clusters["cluster_"+std::to_string(j)] = WriteNode(n);
					j++;
				}
			}
			trees["tree_"+std::to_string(t)] = clusters;
			//reset trees Json object so values aren't carried over to unfilled levels
			//clusters.clear();
		}
		if(_verb > 2) cout << npts << " points at level " << l << endl;
		levels["level_"+std::to_string(l)] = trees;
	}
	_root["levels"] = levels;
	//tree_map.clear();
	return levels;
}


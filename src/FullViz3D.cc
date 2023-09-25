#include "FullViz3D.hh"

FullViz3D::FullViz3D(vector<node*> nodes){
	_nodes = nodes;
	
	for(int i = 0; i < (int)_nodes.size(); i++){
		m_points->AddPoints(*_nodes[i]->points);	
		_k.push_back(_nodes[i]->model->GetNClusters());
	}
	m_n = m_points->GetNPoints();
	

}



json FullViz3D::WriteNode(node* node){
	
	json subcluster = json::object();
	json subclusters = json::object();
	json cluster = json::object();
	json data = json::object();

	vector<double> color;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> w;
	
	vector<double> eigenVec_0;
	vector<double> eigenVec_1;
	vector<double> eigenVec_2;
	
	int kmax = node->model->GetNClusters();
	BasePDFMixture* model = node->model;
	PointCollection* points = node->points;
	if(points->GetNPoints() == 0){
		return cluster;
	}


	//for writing individual cluster to its own plot/its history to its own json
	if(_localcoords){
		//center points on mean (not weighted) such that avg = 0 -> this is a deltaX, deltaY, deltaT plot/coords
	}
	
	for(int i = 0; i < points->GetNPoints(); i++){
		//eta
		x.push_back(points->at(i).Value(0));
		//phi
		y.push_back(points->at(i).Value(1));
		//time
		z.push_back(points->at(i).Value(2));
		//weight
		w.push_back(points->at(i).Weight());

	}
	data["x"] = x;
	data["y"] = y;
	data["z"] = z;
	data["w"] = w;
	cluster["data"] = data;
	

	//set coords for parameter circles
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	vector<double> avgs;
	model->GetAvgWeights(avgs);
	
	double x0, y0, z0;	
	map<string, Matrix> cluster_params;
	for(int k = 0; k < kmax; k++){
		cluster_params = model->GetParameters(k);
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
		subcluster["mixing_coeff_norm"] = cluster_params["pi"].at(0,0);
		subcluster["mu_x"] = x0;
		subcluster["mu_y"] = y0;
		subcluster["mu_z"] = z0;

		subcluster["eigenVal_0"] = eigenVals[0];	
		subcluster["eigenVal_1"] = eigenVals[1];	
		subcluster["eigenVal_2"] = eigenVals[2];	
		
		subcluster["eigenVec_0"] = eigenVec_0;	
		subcluster["eigenVec_1"] = eigenVec_1;	
		subcluster["eigenVec_2"] = eigenVec_2;	
		
		subcluster["color"] = avgs[k];
	
		subclusters["subcluster_"+std::to_string(k)] = subcluster;


		eigenVec_0.clear();
		eigenVec_1.clear();
		eigenVec_2.clear();
	}
	cluster["subclusters"] = subclusters;
	return cluster;
}



void FullViz3D::orderTree(node* node, int level, map<int, NodeStack> &map){
	if(node->val == -1) return;
	map[level].push(node);
	
	orderTree(node->l, level + 1, map);
	orderTree(node->r, level + 1, map);

}

//level{ tree_0{ }, tree_1{ }, ...}
json FullViz3D::WriteLevels(){
	json levels = json::object();
	json trees = json::object();
	json clusters = json::object();
	
	
	int nTrees = (int)_nodes.size(); 

	//create a node - level map for each tree
	vector<map<int, NodeStack>> tree_maps;
	for(int i = 0; i < nTrees; i++){
		map<int, NodeStack> tree_map;
		orderTree(_nodes[i], 0, tree_map);
		tree_maps.push_back(tree_map);	
	}
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
			if(l <= tree_maps[t].rbegin()->first && !tree_maps[t][l].empty()){
				if(_verb > 2) cout << "Tree " << t << ": " << endl;
				node* n = tree_maps[t][l].pop();
				int j = 0;
				//loop through nodes (clusters) at this level for this tree
				while(n->val != -999 && !tree_maps[t][l].empty()){
					if(_verb > 1) cout << "node " << j << " - number of points: " << n->points->GetNPoints() << endl; 
					//if there is a cluster with one point in tree_maps[t][l] (a leaf) add it to tree_maps[t][l+1]
					if(n->points->GetNPoints() == 1 && l <= tree_maps[t].rbegin()->first){
						tree_maps[t][l+1].push(n);
					}
					npts += n->points->GetNPoints();
					clusters["cluster_"+std::to_string(j)] = WriteNode(n);
					n = tree_maps[t][l].pop();
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
	return levels;
}

json FullViz3D::WriteTree(node* root){
	json level = json::object();
	json clusters = json::object();
	map<int, NodeStack> tree_map;
	orderTree(root, 0, tree_map);	

	for (int i = 0; i < tree_map.size(); i++)
		{
		     cout << "Level " << i << ": " << endl;
		     tree_map[i].Print();
		     //for each node in tree_map[i] (NodeStack):
		     node* t = tree_map[i].pop();
		     int j = 0;
		     while(t->val != -999){ 
		     	clusters["cluster_"+std::to_string(j)] = WriteNode(t);
		     	t = tree_map[i].pop();
		     	j++;	
		     }
		     level["clusters_level_"+std::to_string(i)] = clusters;
		}

	json levels = json::object();
	levels["levels"] = level;
	return levels;
}

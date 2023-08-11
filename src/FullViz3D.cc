#include "FullViz3D.hh"

FullViz3D::FullViz3D(vector<node*> nodes){
	_nodes = nodes;
	
	for(int i = 0; i < (int)_nodes.size(); i++){
		m_points->AddPoints(*_nodes[i]->points);	
		_k.push_back(_nodes[i]->model->GetNClusters());
	}
	m_n = m_points->GetNPoints();
	

}



Json::Value FullViz3D::WriteNode(node* node){
	Json::Value cluster;
	Json::Value subclusters;
	Json::Value subcluster;
	Json::Value data;

	Json::Value x(Json::arrayValue);
	Json::Value y(Json::arrayValue);
	Json::Value z(Json::arrayValue);
	Json::Value color;
	
	Json::Value eigenVec_0(Json::arrayValue);
	Json::Value eigenVec_1(Json::arrayValue);
	Json::Value eigenVec_2(Json::arrayValue);
	
	
	int kmax = node->model->GetNClusters();
	BasePDFMixture* model = node->model;
	PointCollection* points = node->points;
	if(points->GetNPoints() == 0){
		return cluster;
	}
	for(int i = 0; i < points->GetNPoints(); i++){
		//eta
		x.append(points->at(i).Value(0));
		//phi
		y.append(points->at(i).Value(1));
		z.append(points->at(i).Value(2));
	}
	data["x"] = x;
	data["y"] = y;
	data["z"] = z;
	data["color"] = node->color;

	cluster["data"] = data;
	//if no points - empty plot
//	if(x.size() == 0) return;

	//set coords for parameter circles
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	
	double x0, y0, z0;	
	map<string, Matrix> cluster_params;
	for(int k = 0; k < kmax; k++){
		cluster_params = model->GetParameters(k);
		x0 = cluster_params["mean"].at(0,0);
		y0 = cluster_params["mean"].at(1,0);
		z0 = cluster_params["mean"].at(2,0);

		cluster_params["cov"].eigenCalc(eigenVals, eigenVecs);
		for(int i = 0; i < 3; i++){
			eigenVec_0.append(eigenVecs[0].at(i,0));
			eigenVec_1.append(eigenVecs[1].at(i,0));
			eigenVec_2.append(eigenVecs[2].at(i,0));
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
	
		subclusters["subcluster_"+std::to_string(k)] = subcluster;


		eigenVec_0.clear();
		eigenVec_1.clear();
		eigenVec_2.clear();
	}
	cluster["subclusters"] = subclusters;
	return cluster;
}






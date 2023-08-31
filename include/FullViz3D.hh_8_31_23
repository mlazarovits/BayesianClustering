#ifndef FullViz3D_HH
#define FullViz3D_HH

#include "ClusterVizBase.hh"
#include "BaseTree.hh"
#include "NodeStack.hh"
#include "json/json.h"
#include <string>
#include <iostream>
#include <fstream>
using std::string;


using node = BaseTree::node;
class FullViz3D : public ClusterVizBase{
	public:
		FullViz3D() : ClusterVizBase(){ };
		FullViz3D(vector<node*> nodes);
		virtual ~FullViz3D(){ };

		//writes individual node
		Json::Value WriteNode(node* node);
		//writes tree following from given root node
		Json::Value WriteTree(node* root);
		void AddPlot(string filename = "test"){ };
		void Write(){ };
		Json::Value WriteLevels();
		void Write(string filename){
			//writing level{ tree{} }
			WriteLevels();	

			Json::StreamWriterBuilder builder;
			const std::string json_file = Json::writeString(builder, _root);

			std::ofstream file;
			file.open(filename+".json");
			file << json_file << endl;
			cout << "Writing to: " << filename << ".json" << endl;
		}
		void SeeData(){ };
	private:
		vector<node*> _nodes;
		vector<int> _k; //number of subclusters per cluster
		//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
		Json::Value _root;

		void orderTree(node* node, int level, map<int, NodeStack> &map);
};
#endif

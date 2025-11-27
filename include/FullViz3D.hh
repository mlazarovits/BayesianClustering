#ifndef FullViz3D_HH
#define FullViz3D_HH

#include "ClusterVizBase.hh"
#include "BaseTree.hh"
#include "NodeStack.hh"
#include "nlohmann/json.hpp"
#include "Jet.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip> 
using std::string;

using json = nlohmann::json;
using node = BaseTree::node;
class FullViz3D : public ClusterVizBase{
	public:
		FullViz3D() : ClusterVizBase(){ };
		virtual ~FullViz3D(){ };

		//writes individual node
		json WriteNode(node *n);
		//writes tree following from given root node
		//json WriteTree(node* root);
		void AddPlot(string filename = "test"){ };
		void Write(){ };
		json WriteLevels(vector<std::shared_ptr<node>>& nodes);
		void Write(vector<std::shared_ptr<node>>& nodes, string filename){
			//writing level{ tree{} }
			WriteLevels(nodes);	
			std::ofstream file;
			file.open(filename+".json");
			////sets 4 space indent
			cout << "Writing to: " << filename << ".json" << endl;
			file << std::setw(4) << _root << endl;
			file.close();
		}
		void SeeData(){ };


		void AddTrueJets(vector<Jet> jets){
			json true_jets = json::object();
			for(int j = 0; j < jets.size(); j++){
				//eta
				true_jets["jet_"+std::to_string(j)]["x"] = jets[j].eta();
				//phi
				true_jets["jet_"+std::to_string(j)]["y"] = jets[j].phi_02pi();
				//time - set right now to be energy weighted time of all rechits (see Jet.hh)
				true_jets["jet_"+std::to_string(j)]["z"] = jets[j].time();
			}
			_root["true_jets"] = true_jets;
	
		}

	private:
		//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
		json _root = json::object();

		void orderTree(node* n, int level, map<int, NodeStack> &map);

};
#endif

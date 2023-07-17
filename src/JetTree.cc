#include "JetTree.hh"
#include <stdexcept>

using std::logic_error;

JetTree::JetTree(){ };

JetTree::~JetTree(){
	m_tree.clear();
}

void JetTree::AddLayer(vector<Jet> jets){
	m_tree.push_back(jets);
}

vector<Jet> JetTree::at(int i){
	return m_tree[i];
}

Jet JetTree::at(int i, int j){
	return m_tree[i][j];
}


int JetTree::GetNLayers(){
	return (int)m_tree.size();
}

int JetTree::GetNJets(){
	int cnt = 0;
	for(int i = 0; i < (int)m_tree.size(); i++)
		cnt += (int)m_tree[i].size();
	return cnt;
}



bool JetTree::check_jet_in_layer(Jet& jet, int i){
	for(int j = 0; j < (int)m_tree[i].size(); j++){
		if(jet == m_tree[i][j])
			return true;
	}
	return false;
}

void JetTree::check_validity(){
	//check that there are less than or equal jets in layer i compared to layer i+1
	Jet p1, p2;
	for(int i = 1; i < (int)m_tree.size(); i++){
		if((int)m_tree[i].size() > (int)m_tree[i-1].size()){
			throw logic_error("Error: layer "+std::to_string(i)+" has more jets than layer "+std::to_string(i-1)+" in JetTree.");
			return;
		}
		
		//check that parents of jet j are in previous layer
		for(int j = 0; j < (int)m_tree[i].size(); j++){
			m_tree[i][j].GetParents(p1, p2);	
			if(!check_jet_in_layer(p1, i-1) || !check_jet_in_layer(p2, i-1)){
				throw logic_error("Error: parent of jet "+std::to_string(j)+" not in layer "+std::to_string(i-1)+" in JetTree.");
				return;
			}
		}


	}

}

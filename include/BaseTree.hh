#ifndef BaseTree_HH
#define BaseTree_HH

#include "PointCollection.hh"

class BaseTree{
	public:
		BaseTree(){ }
		virtual ~BaseTree(){ delete _head; delete _z; delete _t; };
		//node structure
		struct node{
			//points at node or info
			PointCollection* points;
			//left and right subnodes
			struct node *l, *r;
			//posterior value - info or key
			double val;
			//factor in prior
			double d;
			//probability of being in tree T_k p(D_k | T_k)
			double prob_tk;
			//node id? int idx;
		};

		node* _head, *_z, *_t;

};
#endif





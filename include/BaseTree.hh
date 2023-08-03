#ifndef BaseTree_HH
#define BaseTree_HH

#include "PointCollection.hh"

class BaseTree{
	public:
		BaseTree(){
			_z = (struct node*) malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1; _z->d = -1; _z->prob_tk = -1;
			_head = (struct node*)malloc(sizeof *_head);
			_head->r = _z; _head->val = 0;
		}
		virtual ~BaseTree(){ delete _head; delete _z; }
		//node structure
		struct node{
			//points at node or info
			PointCollection* points;
			//left and right subnodes
			struct node *l;
			struct node *r;
			//posterior value - info or key
			double val;
			//factor in prior
			double d;
			//probability of being in tree T_k p(D_k | T_k)
			double prob_tk;
			//for debugging - making sure correct merges are happening	
			//std::string name;
		};

		node* _head, *_z, *_t;

};
#endif





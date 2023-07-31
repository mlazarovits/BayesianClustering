#ifndef SearchTree_HH
#define SearchTree_HH

#include "PointCollection.hh"
#include "BaseTree.hh"
using node = BaseTree::node;

class SearchTree : public BaseTree{
	public:
		SearchTree(){
			_z = (struct node*) malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1.;
			_head = (struct node*) malloc(sizeof *_head);
			//head should have smallest value (therefore nothing in left subtree)
			_head->r = _z; _head->val = 0.; _head->idxs = std::make_pair(-1,-1);
		}
		virtual ~SearchTree(){ };

		//for binary search tree, all records with smaller keys are in the left subtree 
		//and all records in right subtree have larger or equal key values
		void insert(double v, PointCollection* pc){
			struct node *p, *x;
			p = _head; x = _head->r;
			//search while nodes are unterminated
			while(x != _z){
				//check if p is greater than or equal to (right subtree) or less than (left subtree) current node
				p = x; x = (v < x->val) ? x->l : x->r; 
			}
			//once terminated node is found, insert new node
			x = (struct node*) malloc(sizeof *x);
			x->val = v; x->points = pc; x->l = _z; x->r = _z;
			if(v < p->val) p->l = x; else p->r = x;

		}

		void insert(double v, pair<int, int> idxs){
			struct node *p, *x;
			p = _head; x = _head->r;
			//search while nodes are unterminated
			while(x != _z){
				//check if p is greater than or equal to (right subtree) or less than (left subtree) current node
				p = x; x = (v < x->val) ? x->l : x->r; 
			}
			//once terminated node is found, insert new node
			x = (struct node*) malloc(sizeof *x);
			x->val = v; x->idxs = idxs; x->l = _z; x->r = _z;
			if(v < p->val) p->l = x; else p->r = x;

		}

	
		//search for given posterior value, return idxs of nodes that make up merge
		pair<int, int> search(double v){
			struct node *x = _head->r;
			//stops search to end in termination node
			_z->val = v;
			while(v != x->val){
				x = (v < x->val) ? x->l : x->r;
			}
			return x->idxs;
		}	
	
		//delete based on node idxs	
		void del(pair<int, int> idxs){
			struct node *c, *p, *x;
			_z->idxs = idxs;
			p = _head; x = _head->r;
			double v = x->val;
			while(idxs.first != x->idxs.first && idxs.second != x->idxs.second){
				p = x;
				x = (v < x->val) ? x->l : x->r;
			}
			_t = x;
			v = _t->val;
			if(_t->r == _z) x = x->l;
			else if(_t->r->l == _z){
				x = x->r;
				x->l = _t->l;
			}
			else{
				c = x->r;
				while(c->l->l != _z) c = c->l;
				x = c->l; c->l = x->r;
				x->l = _t->l; x->r = _t->r;
			}
			free(_t);
			if(v < p->val) p->l = x;
			else p->r = x;
		}

		//delete based on node idxs	
		void del(int idx){
			struct node *c, *p, *x;
			_z->idxs.first = idx;
			p = _head; x = _head->r;
			double v = x->val;
			while(idx != x->idxs.first || idx != x->idxs.second){
				p = x;
				x = (v < x->val) ? x->l : x->r;
			}
			_t = x;
			v = _t->val;
			if(_t->r == _z) x = x->l;
			else if(_t->r->l == _z){
				x = x->r;
				x->l = _t->l;
			}
			else{
				c = x->r;
				while(c->l->l != _z) c = c->l;
				x = c->l; c->l = x->r;
				x->l = _t->l; x->r = _t->r;
			}
			free(_t);
			if(v < p->val) p->l = x;
			else p->r = x;
		}

		void del(double v){
			struct node *c, *p, *x;
			_z->val = v;
			p = _head; x = _head->r;
			while(v != x->val){
				p = x;
				x = (v < x->val) ? x->l : x->r;
			}
			_t = x;
			if(_t->r == _z) x = x->l;
			else if(_t->r->l == _z){
				x = x->r;
				x->l = _t->l;
			}
			else{
				c = x->r;
				while(c->l->l != _z) c = c->l;
				x = c->l; c->l = x->r;
				x->l = _t->l; x->r = _t->r;
			}
			free(_t);
			if(v < p->val) p->l = x;
			else p->r = x;
		}
 



};
#endif

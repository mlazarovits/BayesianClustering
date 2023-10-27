#ifndef BaseTree_HH
#define BaseTree_HH

#include "PointCollection.hh"
#include "BasePDFMixture.hh" 

class BaseTree{
	public:
		BaseTree(){
			_z = (struct node*) malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1; _z->d = -1; _z->prob_tk = -1; _z->model = nullptr; //_z->color = -1; 
			_z->points = new PointCollection();
			_z->idx = -1;

			_head = (struct node*)malloc(sizeof *_head);
			_head->l = _head; _head->r = _head; _head->val = -1; _head->d = -1; _head->prob_tk = -1; _head->model = nullptr;
			_head->points = new PointCollection();
			_head->idx = -1;
		}
		virtual ~BaseTree(){ free(_head); free(_z); }
		//node structure
		struct node{
			//points at node or info
			PointCollection* points = nullptr;
			//left and right subnodes
			struct node *l = nullptr;
			struct node *r = nullptr;
			//posterior value - info or key
			double val;
			//factor in prior
			double d;
			//model of cluster
			BasePDFMixture* model = nullptr;
			//probability of being in tree T_k p(D_k | T_k)
			double prob_tk;
			//3D distance to nearest neighbor
			double nn3dist;
			//mirror node
			struct node* mirror = nullptr;
			int idx; //index of node in merge tree
			//color for plotting -> maps to color map
			//int color;
			//for debugging - making sure correct merges are happening	
			//std::string name;
			
			//node* operator =(const node* n){
			node(const node& n){
				points = n.points;
				l = n.l;
				r = n.r;
				val = n.val;
				d = n.d;
				model = n.model;
				prob_tk = n.prob_tk;
			}

			//destructor for node
			~node(){
				delete points;
				delete model;
				free(l);
				free(r);
				free(mirror);
			}

		};
		struct listnode{
			//posterior value in here
			node* n;
			struct listnode* next;
		};


		node* _head, *_z, *_t;
};
#endif





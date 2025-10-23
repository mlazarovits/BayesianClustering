#ifndef BaseTree_HH
#define BaseTree_HH

#include "PointCollection.hh"
#include "GaussianMixture.hh" 
//#include <boost/multiprecision/cpp_bin_float.hpp>
//using namespace boost::multiprecision;

class BaseTree{
	public:
		BaseTree(){
			_z = (struct node*) malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1; 
			//_z->d = -1; _z->prob_tk = -1; 
			_z->log_val = -1; _z->log_d = -1; _z->log_prob_tk = -1; 
			_z->model = nullptr;
			_z->points = new PointCollection();
			_z->idx = -1;
			_z->nndist = -999;

			_head = (struct node*)malloc(sizeof *_head);
			_head->l = _head; _head->r = _head; _head->val = -1; 
			//_head->d = -1; _head->prob_tk = -1; 
			_head->log_val = -1; _head->log_d = -1; _head->log_prob_tk = -1; 
			_head->model = nullptr;
			_head->points = new PointCollection();
			_head->idx = -1;
			_head->nndist = -999;
		}
		virtual ~BaseTree(){ free(_head); free(_z); }
		//node structure
		struct node{
			//points at node or info
			PointCollection* points = new PointCollection();//nullptr;
			//left and right subnodes
			struct node *l = nullptr;
			struct node *r = nullptr;
			//posterior value - info or key
			double val = -999;
			double log_val = -999;
			//log(p_dk_h1*pi)
			double log_h1_prior = -1e100;
			//log(p_di_ti*p_dj_tj*di*dj)
			double log_didj = -1e100;
			//factor in prior
			//cpp_bin_float_100 d = -999;
			double log_d = -999;
			//model of cluster
			GaussianMixture* model = nullptr;
			//probability of being in tree T_k p(D_k | T_k)
			//double prob_tk = -999;
			//cpp_bin_float_100 prob_tk = -999;
			double log_prob_tk = -999;
			//3D distance to nearest neighbor
			double nndist = 1e300;
			//mirror node
			struct node* mirror = nullptr;
			//is mirror
			bool ismirror = false;
			int idx = -999; //index of node in merge tree
			
			//node* operator =(const node* n){
			node(const node& n){
			//	if(n.points)
			//		points = std::make_unique<PointCollection>(*n.points);
			//	else
			//		points = nullptr;
				points = n.points;
				l = n.l;
				r = n.r;
				val = n.val;
				log_val = n.log_val;
				log_h1_prior = n.log_h1_prior;
				log_didj = n.log_didj;
				//d = n.d;
				log_d = n.log_d;
				model = n.model;
				//prob_tk = n.prob_tk;
				log_prob_tk = n.log_prob_tk;
				nndist = n.nndist;
				ismirror = n.ismirror;
				mirror = n.mirror;
				idx = n.idx;
			}

			//destructor for node
			~node(){
				//cout << "~node start" << endl;
				//cout << "b" << endl;
				model = nullptr;
				//cout << "c" << endl;
				delete model;
				//cout << "d" << endl;
				//don't free l + r because those could be freed themselves as nodes
				if(mirror != NULL) free(mirror);
				//cout << "~node done" << endl;
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





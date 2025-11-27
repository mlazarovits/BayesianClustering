#ifndef BaseTree_HH
#define BaseTree_HH

#include "PointCollection.hh"
#include "GaussianMixture.hh" 

class BaseTree{
	public:
		BaseTree(){ }
		virtual ~BaseTree(){  }
		//node structure
		struct node{
			//points at node or info
			std::unique_ptr<PointCollection> points = std::make_unique<PointCollection>();//nullptr;
			//left and right subnodes
			//std::shared_ptr<node> l = nullptr;
			//std::shared_ptr<node> r = nullptr;
			node* l = nullptr;
			node* r = nullptr;
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
			std::unique_ptr<GaussianMixture> model = nullptr;
			//probability of being in tree T_k p(D_k | T_k)
			//double prob_tk = -999;
			//cpp_bin_float_100 prob_tk = -999;
			double log_prob_tk = -999;
			//3D distance to nearest neighbor
			double nndist = 1e300;
			//mirror node
			node* mirror = nullptr;
			//is mirror
			bool ismirror = false;
			int idx = -999; //index of node in merge tree
			
			node(){ };

			node(const node& n):
				val(n.val),
				log_val(n.log_val),
				log_h1_prior(n.log_h1_prior),
				log_didj(n.log_didj),
				log_d(n.log_d),
				log_prob_tk(n.log_prob_tk),
				nndist(n.nndist),
				ismirror(n.ismirror),
				idx(n.idx),
				//pointers
				points(n.points ? std::make_unique<PointCollection>(*n.points) : nullptr),	
				model(n.model ? std::make_unique<GaussianMixture>(*n.model) : nullptr),	
				l(n.l),	
				r(n.r),	
				mirror(n.mirror)	
			{ }


			node& operator=(const node& n){
				//can't self assign
				if(this != &n){
					val = n.val;
					log_val = n.log_val;
					log_h1_prior = n.log_h1_prior;
					log_didj = n.log_didj;
					log_d = n.log_d;
					log_prob_tk = n.log_prob_tk;
					nndist = n.nndist;
					ismirror = n.ismirror;
					idx = n.idx;
				
					points = n.points ? std::make_unique<PointCollection>(*n.points) : nullptr;	
					model = n.model ? std::make_unique<GaussianMixture>(*n.model) : nullptr;	
					mirror = n.mirror;
					l = n.l;	
					r = n.r;	
				}
				return *this;
			}

			//destructor for node
			~node(){ }

		};

};
#endif





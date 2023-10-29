#ifndef NodeStack_HH
#define NodeStack_HH

#include "RandomSample.hh"
#include "PointCollection.hh"
#include "BaseTree.hh"

using node = BaseTree::node;
using listnode = BaseTree::listnode;
class NodeStack{
	public:
		NodeStack(){
			_head = (listnode*)malloc(sizeof *_head);
			node* h = (node*)malloc(sizeof* h); h->val = 999; //sort high to low
			h->l = h; h->r = h; h->d = -1; h->prob_tk = -1; h->model = nullptr; //h->color = -1; 
			h->points = new PointCollection();
			_head->n = h;	
		
			_z = (listnode*)malloc(sizeof *_z);
			node* z = (node*)malloc(sizeof* z); z->val = -999; //sort high to low
			z->l = z; z->r = z; z->d = -1; z->prob_tk = -1; z->model = nullptr; //z->color = -1; 
			z->points = new PointCollection(); 
			_z->n = z;
			
			_head->next = _z; _z->next = _z; 
		}		

		NodeStack(const NodeStack& nodes){
			_head = nodes._head;
			_z = nodes._z;
		}
		virtual ~NodeStack(){ };//free(_head); free(_z); };


		listnode* GetList(){
			return _head;
		}

		//inserts x after t	
		void insertafter(listnode* t, listnode* x){
			x->next = t->next; t->next = x;
		}

		void insert(node* x){
			listnode* y = (listnode*)malloc(sizeof *y);
			y->n = x;
			insertafter(_head,y);	
		}

		void push(node* x){
			listnode* t = (listnode *)malloc(sizeof *t);
			t->n = x; t->next = _head->next;
			_head->next = t;
		}

		//deletes node after t
		void deletenext(listnode *t){
			t->next = t->next->next;
		}

		//get node from top of stack (should be largest)
		node* pop(){
			if(empty()) return nullptr;
			node *x = (node*)malloc(sizeof *x);
			listnode* t = _head->next; _head->next = t->next;
			x = t->n; 
			free(t);
			return x;	
		}

		//pops off top node and removes then impossible nodes (merges with one of the parents of given merge)
		node* fullpop(){
			if(empty()) return nullptr;
			node *x = pop();
			listnode *c;
			c = _head;
			while(c->next != _z){
				//cout << "looking at node: " << c->next->node->val << endl;
				//remove nodes whose parents (either l or r) are in the max merge
				if(c->next->n->l == x->l || c->next->n->r == x->r || c->next->n->l == x->r || c->next->n->r == x->l){
					deletenext(c);	
				}
				//update c to next listnode
				else c = c->next;
			}
			return x;	
		}
		

		void merge(const NodeStack& list){
			listnode* a = list._head->next;
			_head = merge(_head, a);
		}

 

		listnode* merge(listnode* a, listnode* b){
			listnode* c;
			c = _z;
			do{
				if(a->n->val >= b->n->val){
      				c->next = a; c = a; a = a->next;
      			}
      			else{
      				c->next = b; c = b; b = b->next;
      			}
      			}while(c != _z);
      			c = _z->next; _z->next = _z;
			return c;
		}


		listnode* mergesort(listnode* c){
			//cout << "NodeStack::mergesort" << endl;
			//cout << "c post: " << c->node->val << endl;
			listnode* a, *b;
			if(c->next != _z){
				a = c; b = c->next->next->next;
				while(b != _z){
					c = c->next; b = b->next->next;
				}
				b = c->next; c->next = _z;
				//cout << "NodeStack::mergesort - end 1" << endl;
		      	return merge(mergesort(a), mergesort(b));
		      }	
			//cout << "NodeStack::mergesort - end 2 - c->next: " << c->next->node->val << endl;
			return c;
		}
		

		void sort(){
		//listnode* c = _head->next;
		//cout << "unsorted" << endl;	
		//int i = 0;
		//while(c != _z){ cout << i << " " << c->n->val << endl; i++; c = c->next; } 
		_head = mergesort(_head);
		//cout << "sorted" << endl;	
		//listnode* g = _head;
		//i = 0;
		//while(g != _z){ cout << i << " " << g->n->val << endl; i++; g = g->next; } 
		}

		bool empty(){
			return _head->next == _z; 
		}

		void Print(int v = 0){
			if(empty()) return;
			listnode* g = _head->next;
			int i = 1;
			while(g != _z){ 
				if(v == 1){
					cout << "nodes " << g->n->l->idx << ": ";
					g->n->l->points->Print();
					cout << " and " << g->n->r->idx << ": ";
					g->n->r->points->Print();
				}
				if(v == 0){
					cout << "cluster " << i << ": ";
				} 
				cout << "this rk: " << g->n->val << "\n" << endl;
			i++; g = g->next; } 
		}

		private:
			listnode* _head = nullptr; 
			listnode* _z = nullptr;



};
#endif

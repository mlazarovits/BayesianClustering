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
			_z = (listnode*)malloc(sizeof *_z);
			_head->next = _z; _z->next = _z; 
			node* h = (node*)malloc(sizeof* h); h->val = 999; //sort high to low
			_head->node = h;
			node* z = (node*)malloc(sizeof* z); z->val = -999; //sort high to low
			_z->node = z; 
		}		

		NodeStack(const NodeStack& nodes){
			_head = nodes._head;
			_z = nodes._z;
		}
		virtual ~NodeStack(){ };


		listnode* GetList(){
			return _head;
		}

		//inserts x after t	
		void insertafter(listnode* t, listnode* x){
			x->next = t->next; t->next = x;
		}

		void insert(node* x){
			listnode* y = (listnode*)malloc(sizeof *y);
			y->node = x;
			insertafter(_head,y);	
		}

		void push(node* x){
			listnode* t = (listnode *)malloc(sizeof *t);
			t->node = x; t->next = _head->next;
			_head->next = t;
		}

		//deletes node after t
		void deletenext(listnode *t){
			t->next = t->next->next;
		}

		//get node from top of stack (should be largest)
		node* pop(){
			node *x;
			listnode *t;
			t = _head->next; _head->next = t->next;
			x = t->node;
			free(t);
			return x;	
		}

		//pops off top node and removes then impossible nodes (merges with one of the parents of given merge)
		node* fullpop(){
			node *x = pop();
			listnode *c;
			c = _head;
			while(c->next != _z){
				//cout << "looking at node: " << c->next->node->val << endl;
				//remove nodes whose parents (either l or r) are in the max merge
				if(c->next->node->l == x->l || c->next->node->r == x->r || c->next->node->l == x->r || c->next->node->r == x->l){
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
			//cout << "NodeStack::merge 1" << endl;
			c = _z;
			//cout << "NodeStack::merge 2 " << a->node->val << endl;
			do{
			//	cout << "begin do-while loop - a: " << a->node->val << " b: " << b->node->val << endl;
				if(a->node->val >= b->node->val){
			//cout << "NodeStack::merge 3" << endl;
      				c->next = a; c = a; a = a->next;
			//cout << "NodeStack::merge 4" << endl;
      			}
      			else{
			//cout << "NodeStack::merge 5" << endl;
      				c->next = b; c = b; b = b->next;
			//cout << "NodeStack::merge 6" << endl;
      			}
			//cout << "NodeStack::merge 7" << endl;
      		}while(c != _z);
			//cout << "NodeStack::merge 8" << endl;
      		c = _z->next; _z->next = _z;
			//cout << "NodeStack::merge - end" << endl;
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
		//struct listnode* c = _head->next;
		//cout << "unsorted" << endl;	
		//int i = 0;
		//while(c != _z){ cout << i << " " << c->node->val << endl; i++; c = c->next; } 
		_head = mergesort(_head);
		//cout << "sorted" << endl;	
		//struct listnode* g = _head;
		//i = 0;
		//while(g != _z){ cout << i << " " << g->node->val << endl; i++; g = g->next; } 
		}


		void Print(){
			listnode* g = _head->next;
			int i = 1;
			while(g != _z){ cout << i << " " << g->node->val << endl; i++; g = g->next; } 
		}

		private:
			listnode* _head, *_z;



};
#endif

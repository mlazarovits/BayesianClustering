#ifndef NodeList_HH
#define NodeList_HH

#include "RandomSample.hh"
#include "PointCollection.hh"
#include "BaseTree.hh"

using node = BaseTree::node;
class NodeList{
	public:
		NodeList(){
			_head = (struct listnode*)malloc(sizeof *_head);
			_z = (struct listnode*)malloc(sizeof *_z);
			_head->next = _z; _z->next = _z;
		}		

		NodeList(const NodeList& nodes){
			_head = nodes._head;
			_z = nodes._z;
		}
		virtual ~NodeList(){ };
		struct listnode{
			//posterior value in here
			node* node;
			struct listnode* next;
		};


		listnode* GetList(){
			return _head;
		}

		//inserts x after t	
		void insertafter(listnode* t, listnode* x){
			x = (struct listnode*) malloc(sizeof *x);
			x->next = t->next; t->next = x;
		}

		void insert(listnode* x){
			insertafter(_head, x); 
		}

		void insert(node* x){
			struct listnode* y = (struct listnode*)malloc(sizeof *y);
			y->node = x;
			insert(y);	
		}

		void push(node* x){
			struct listnode* t = (struct listnode *)malloc(sizeof *t);
			t->node = x; t->next = _head->next;
			_head->next = t;
		}

		//deletes node after t
		void deletenext(struct listnode *t){
			t->next = t->next->next;
		}

		//get node from top of stack (should be largest)
		node* pop(){
			node *x;
			struct listnode *t;
			t = _head->next; _head->next = t->next;
			x = t->node;
			free(t);
			return x;	
		}

		//pops off top node and removes then impossible nodes (merges with one of the parents of given merge)
		node* fullpop(){
			node *x;
			struct listnode *t, *c;
			t = _head->next; _head->next = t->next;
			x = t->node;
			c = _head->next;
			while(c != _z){
				if(c->next->node->l ==  x->l || c->next->node->r == x->r){
					deletenext(c);	
				}
			}
			free(t);
			return x;	

		}


		//TODO: check
		void merge(const NodeList& list1, const NodeList& list2){
			struct listnode* a = list1._head;
			struct listnode* b = list2._head;
			struct listnode* c = _head;	
			if(c->next != _z){
				a = c; b = c->next->next->next;
				while(b != _z){
					c = c->next; b = b->next->next;
				}
				b = c->next; c->next = _z;
				_head = merge(mergesort(a), mergesort(b));
			}	
		}
 

		struct listnode* merge(struct listnode* a, struct listnode* b){
			struct listnode* c;
			c = _z;
			do
				if(a->node->val <= b->node->val){
					c->next = a; c = a; a = a->next;
				}
				else{
					c->next = b; c = b; b = b->next;
				}
			while(c != _z);
			c = _z->next; _z->next = _z;
			return c;
		}


		struct listnode* mergesort(listnode* c){
			struct listnode* a, *b;
			if(c->next != _z){
				a = c; b = c->next->next->next;
				while(b != _z){
					c = c->next; b = b->next->next;
				}
				b = c->next; c->next = _z;
				return merge(mergesort(a), mergesort(b));
			}	
			return c;
		}
		

		void sort(){
			_head = mergesort(_head);
		}

		private:
			listnode* _head, *_z;



};
#endif

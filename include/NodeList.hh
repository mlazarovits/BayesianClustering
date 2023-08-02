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
			node* h = (node*)malloc(sizeof* h); h->val = 999; //sort high to low
			_head->node = h;
			node* z = (node*)malloc(sizeof* z); z->val = -999; //sort high to low
			_z->node = z; 
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
			x->next = t->next; t->next = x;
		}

		void insert(node* x){
			struct listnode* y = (struct listnode*)malloc(sizeof *y);
			y->node = x;
			insertafter(_head,y);	
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
			node *x = pop();
			struct listnode *c;
			c = _head->next;
			while(c != _z){
				if(c->next->node->l == x->l || c->next->node->r == x->r){
					deletenext(c);	
				}
				//update c to next listnode
				else c = c->next;
			}
			return x;	
		}
		

		//TODO: check
		void merge(const NodeList& list){
			struct listnode* a = list._head;
//			struct listnode* c = _head;	
			_head = merge(_head, a);
		/*
			if(c->next != _z){
				a = c; b = c->next->next->next;
				while(b != _z){
					c = c->next; b = b->next->next;
				}
				b = c->next; c->next = _z;
				_head = merge(mergesort(a), mergesort(b));
			}	
		*/
		}

 

		struct listnode* merge(struct listnode* a, struct listnode* b){
			struct listnode* c;
			//cout << "NodeList::merge 1" << endl;
			c = _z;
			//cout << "NodeList::merge 2 " << a->node->val << endl;
			do{
			//	cout << "begin do-while loop - a: " << a->node->val << " b: " << b->node->val << endl;
				if(a->node->val >= b->node->val){
			//cout << "NodeList::merge 3" << endl;
      				c->next = a; c = a; a = a->next;
			//cout << "NodeList::merge 4" << endl;
      			}
      			else{
			//cout << "NodeList::merge 5" << endl;
      				c->next = b; c = b; b = b->next;
			//cout << "NodeList::merge 6" << endl;
      			}
			//cout << "NodeList::merge 7" << endl;
      		}while(c != _z);
			//cout << "NodeList::merge 8" << endl;
      		c = _z->next; _z->next = _z;
			//cout << "NodeList::merge - end" << endl;
			return c;
		}


		struct listnode* mergesort(listnode* c){
			//cout << "NodeList::mergesort" << endl;
			//cout << "c post: " << c->node->val << endl;
			struct listnode* a, *b;
			if(c->next != _z){
				a = c; b = c->next->next->next;
				while(b != _z){
					c = c->next; b = b->next->next;
				}
				b = c->next; c->next = _z;
				//cout << "NodeList::mergesort - end 1" << endl;
		      	return merge(mergesort(a), mergesort(b));
		      }	
			//cout << "NodeList::mergesort - end 2 - c->next: " << c->next->node->val << endl;
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
		struct listnode* g = _head->next;
		int i = 1;
		while(g != _z){ cout << i << " " << g->node->val << endl; i++; g = g->next; } 
		}

		private:
			listnode* _head, *_z;



};
#endif

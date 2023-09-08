#ifndef NodeQueue_HH
#define NodeQueue_HH
#include "BaseTree.hh"

using node = BaseTree::node;
using listnode = BaseTree::listnode;
class NodeQueue{
	public:
		NodeQueue(){
			_head = (listnode*)malloc(sizeof *_head);
			_tail = (listnode*)malloc(sizeof *_tail);
			_none = (listnode*)malloc(sizeof *_none);
		
			//end delimiter
			_z = (node*)malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1; _z->d = -1; _z->prob_tk = -1; _z->model = nullptr; _z->color = -1;
			_none->n = _z;
			_none->next = _none;	

			_head = _tail = _none;
					
		};
		virtual ~NodeQueue(){ };
		NodeQueue(const NodeQueue& nodes){
			_head = nodes._head;
			_tail = nodes._tail;
		};


		//put element at top 
		void put(node* node){
			listnode* t = (listnode*)malloc(sizeof *t);
			t->n = node;
			//if empty
			if(empty()){
				t->next = _none;
				_head = _tail = t; 
				return;
			}
			//put new node at end
			t->next = _none;
			_tail->next = t;	
			_tail = t;
		};
	
		//get _head
		node* get(){
			if(empty()) return nullptr;
			listnode* t = _head;
			_head = _head->next;
			if(_head == _none) _tail = _none;
			node* x = t->n;
			free(t);
			return x;
		
		};


		bool empty(){ return _head == _none; }


		void Print(){
			listnode* t = (listnode*)malloc(sizeof *t);
			t = _head;
			while(t != _none){
				cout << t->n->val << endl;
				t = t->next;
			}

		}

	private:
		listnode* _head, *_tail, *_none;
		node* _z;
};
#endif

#include "BaseTree.hh"

using node = BaseTree::node;
using listnode = BaseTree::listnode;
class NodeQueue{
	public:
		NodeQueue(){
			_head = (listnode*)malloc(sizeof *_head);
			_tail = (listnode*)malloc(sizeof *_tail);
			
			//end delimiter
			_z = (node*)malloc(sizeof *_z);
			_z->l = _z; _z->r = _z; _z->val = -1; _z->d = -1; _z->prob_tk = -1; _z->model = nullptr; _z->color = -1;
	
			_head->node = _z;
			_tail->node = _z;
					
		};
		virtual ~NodeQueue(){ };
		NodeQueue(const NodeQueue& nodes){
			_head = nodes._head;
			_tail = nodes._tail;
		};


		//put element at top 
		void put(node* node){
			t = (listnode*)malloc(sizeof *t);
			t->node = node; t->next = _head->next;
			_head->next = t;
			
		};
	
		//get element behind _tail	
		node* get(){
			listnode* t = _head->next;
			_head->next = t->next;
			node* x = t->node;
			free(t);
			return x;
		
		};


		bool empty(){ return _head == _tail; }

	private:
		listnode* _head, _tail;
		node* _z;
};
#endif

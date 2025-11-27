#ifndef NodeStack_HH
#define NodeStack_HH

#include "RandomSample.hh"
#include "PointCollection.hh"
#include "BaseTree.hh"
#include <memory>

using node = BaseTree::node;
class NodeStack{
	public:
		struct listnode{
			//posterior value in here
			std::shared_ptr<node> n;
			std::unique_ptr<listnode> next;

			listnode(std::shared_ptr<node> in) :
				n(std::move(in)), next(nullptr) { }

		};
		NodeStack() : 
			_head(nullptr)
		{ }		

		//no two nodestacks own the same nodes - breaks single ownership
		NodeStack(const NodeStack& nodes) = delete;
		/*
		{
			_head = std::make_unique<listnode>(*nodes._head);
			//_z = nodes._z;
		}
		*/
		NodeStack& operator=(const NodeStack&) = delete;
		NodeStack(NodeStack&& other) noexcept :
			_head(std::move(other._head)) 
		{ }


		NodeStack& operator=(NodeStack&& other) noexcept {
			if (this != &other) {
				_head = std::move(other._head);
			}
			return *this;
		}


		virtual ~NodeStack(){ };//free(_head); free(_z); };


		listnode* GetList(){
			return _head.get();
		}

		//inserts x after t	
		//void insertafter(listnode* t, listnode* x){
		void insertafter(listnode* t, std::shared_ptr<node> x){
			if(!t) return;
			auto newnode = std::make_unique<listnode>(std::move(x));
			newnode->next = std::move(t->next);
			t->next = std::move(newnode);
		}

		void insert(std::shared_ptr<node> x){
			insertafter(_head.get(),std::move(x));	
		}

		void push(std::shared_ptr<node> x){
			auto newnode = std::make_unique<listnode>(std::move(x));
			newnode->next = std::move(_head);
			_head = std::move(newnode);
		}

		//deletes node after t
		void deletenext(listnode *t){
			if(!t || !t->next) return;
			t->next = std::move(t->next->next);
		}

		//get node from top of stack (should be largest)
		std::shared_ptr<node> pop(){
			if(!_head) return nullptr;
			std::shared_ptr<node> x = std::move(_head->n);
			_head = std::move(_head->next);
			return x;	
		}

		//pops off top node and removes then impossible nodes (merges with one of the parents of given merge)
		std::shared_ptr<node> fullpop(){
			if(empty()) return nullptr;
			auto x = pop();
			auto c = std::move(_head);
			while(c->next != nullptr){
				//cout << "looking at node: " << c->next->node->val << endl;
				//remove nodes whose parents (either l or r) are in the max merge
				if(c->next->n->l == x->l || c->next->n->r == x->r || c->next->n->l == x->r || c->next->n->r == x->l){
					deletenext(c.get());	
				}
				//update c to next listnode
				else c = std::move(c->next);
			}
			return x;	
		}
	

		/*	
		void merge(const NodeStack& list){
			if(!list.head) return;
			listnode* tail = list.head.get();
			while(tail->next){
				tail = tail->next.get();
			}
			tail->next = std::move(head);
			head = std::move(list.head);
		}

 

		unique_ptr<listnode> merge(unique_ptr<listnode> a, unique_ptr<listnode> b){
			auto c = make_unique<listnode>();
			do{
				if(a->n->val >= b->n->val){
      				c->next = a; c = a; a = a->next;
      			}
      			else{
      				c->next = b; c = b; b = b->next;
      			}
      			}while(c != nullptr);
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
		*/

		std::unique_ptr<listnode> merge_sorted(std::unique_ptr<listnode> a, std::unique_ptr<listnode> b) {
		    if (!a) return b;
		    if (!b) return a;
		
		    if (a->n->val <= b->n->val) {
		        a->next = merge_sorted(std::move(a->next), std::move(b));
		        return a;
		    } else {
		        b->next = merge_sorted(std::move(a), std::move(b->next));
		        return b;
		    }
		}
	
		void split_list(std::unique_ptr<listnode> source,
		                std::unique_ptr<listnode>& front,
		                std::unique_ptr<listnode>& back) {
		    if (!source || !source->next) {
		        front = std::move(source);
		        back = nullptr;
		        return;
		    }
		
		    listnode* slow = source.get();
		    listnode* fast = source->next.get();
		
		    while (fast && fast->next) {
		        slow = slow->next.get();
		        fast = fast->next->next.get();
		    }
		
		    front = std::move(source);
		    back = std::move(slow->next);  // cut the list
		}


		std::unique_ptr<listnode> merge_sort_recursive(std::unique_ptr<listnode> node) {
		    if (!node || !node->next) return node; // base case
		
		    std::unique_ptr<listnode> a, b;
		    split_list(std::move(node), a, b);
		
		    a = merge_sort_recursive(std::move(a));
		    b = merge_sort_recursive(std::move(b));
		
		    return merge_sorted(std::move(a), std::move(b));
		}


		void sort(){
			_head = merge_sort_recursive(std::move(_head));
		}
		bool empty(){
			return _head == nullptr; 
		}

		void clear(){
			_head.reset();
		}

		void Print(int v = 0){
			if(empty()) return;
			listnode* g = _head->next.get();
			int i = 1;
			while(g != nullptr){ 
				if(v == 1){
					cout << "nodes " << g->n->l->idx << ": ";
					if(v > 1) g->n->l->points->Print();
					cout << " and " << g->n->r->idx << ": ";
					if(v > 1) g->n->r->points->Print();
				}
				if(v == 0){
					cout << "cluster " << i << ": ";
				} 
				cout << "this rk: " << g->n->val << "\n" << endl;
			i++; g = g->next.get(); } 
		}

		private:
			std::unique_ptr<listnode> _head; 



};
#endif

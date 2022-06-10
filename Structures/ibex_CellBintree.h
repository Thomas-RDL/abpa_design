//============================================================================
//                                  I B E X
// File        : ibexCellBdfs.h
// Author      : Thomas Richard de Latour
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : June, 00
// Last Update : June, 00
//============================================================================
#ifndef __IBEX_BINTREE_H__
#define __IBEX_BINTREE_H__

#include "ibex.h"
#include "Logger.h" // Logger to store and print data in csv

namespace ibex {

struct TNode {
	Cell *cell = NULL;
	double crit = INT_MIN;
	bool test = true;
	inline bool operator==(TNode node){
		if(node.crit == crit && node.cell == cell)
			return true;
		else
			return false;
	}

} ;

bool operator<(const TNode a, const TNode b){
	return a.crit < b.crit;
}

/** \ingroup strategy
 *
 * \brief Cell sorted tree. *
 */
class CellBintree : public CellBuffer {
public:


	CellBintree();

	/** Copy from CellList **/
	/** Flush the buffer.
	 * All the remaining cells will be *deleted* */
	void flush();

	/** Flush the buffer.
	 * All the remaining cells will be *deleted* */
	void flush_tree(TNode* node);

	/** Flush the buffer.
	 * All the remaining cells will be *deleted* */
	void flush_list();

	/** Return the size of the buffer. */
	unsigned int size() const;

	/** Return the size of the tree buffer. */
	unsigned int size_tree() const;

	/** Return true if the buffer is empty. */
	bool empty() const;

	/** push a new cell on the back of the list. */
	void push(Cell* cell);

	/** push a new cell on the front of the list. */
	void push_back(Cell* cell);

	/** Pop a cell from the front of list and return it.*/
	Cell* pop();

	/** Return the next box (but does not pop it).*/
	Cell* top() const;

	/** empty the list in the tree*/
	void insert(const Logger& manif,double (*critere)(const Logger&, CellBuffer&, Cell&));

	/** empty the list in the tree and update criteria for dynamic strategies*/
	void update(const Logger& manif,double (*critere)(const Logger&, CellBuffer&, Cell&));

	/** verify if the list is empty. if YES put the next box from the set in the list*/
	void init_list(bool max=true);

	/** Dynamic behaviour or static behaviour **/
	bool dynamic;

protected :

	/* DFS list Stack of pair of cells  */
	std::list<Cell*> Clist;

	/* BestFS multiset of node of cells and double */
	std::multiset<TNode> Nset;

};

inline CellBintree::CellBintree(){
	dynamic = false;
}

inline 	void CellBintree::insert(const Logger& manif,double (*critere)(const Logger&, CellBuffer&, Cell&)){
	TNode node;
	for(Cell* cell : Clist){
		double crit = (*critere)(manif, *this, *cell);
		node.cell = cell;
		node.crit = crit;
		Nset.insert(node);
	}
	Clist.clear();
	return;
}

inline 	void CellBintree::update(const Logger& manif,double (*critere)(const Logger&, CellBuffer&, Cell&)){
	std::multiset<TNode> _Nset;
	TNode node;

	if(dynamic && Nset.size()>0){
		for(auto it=Nset.begin(); it != Nset.end();it++){
			Logger temp(manif[0].size());
			temp.add(manif[manif.size()-1]);
			temp.set_box(*manif.init_box);
			double crit = (*critere)(temp, *this,(*it->cell));
			node.cell = it->cell;
			if(crit < it->crit)
				node.crit = crit;
			else
				node.crit = it->crit;
			_Nset.insert(node);
		}
		Nset.clear();
		Nset =_Nset;
	}
	for(Cell* cell : Clist){
		double crit = (*critere)(manif, *this, *cell);
		node.cell = cell;
		node.crit = crit;
		Nset.insert(node);
	}
	Clist.clear();
	return;
}

inline void CellBintree::init_list(bool max){
	if(Clist.empty() and !Nset.empty()){
		std::multiset<TNode>::iterator it = Nset.end();
		--it;
		push(it->cell);
		Nset.erase(it);
	}
}


inline void CellBintree::push_back(Cell* cell) {
	if (capacity>0 && size()==capacity) throw CellBufferOverflow();
	Clist.push_back(cell);
}

/** Copy from CelList**/

inline void CellBintree::flush() {
	Nset.clear();
	flush_list();
}

inline void CellBintree::flush_list() {
	while (!Clist.empty()) {
		delete Clist.front();
		Clist.pop_front();
	}
}

inline unsigned int CellBintree::size() const {
	return Clist.size();
}

inline unsigned int CellBintree::size_tree() const{
	return Nset.size();
}

inline bool CellBintree::empty() const {
	return (Clist.empty() && Nset.empty());
}

inline void CellBintree::push(Cell* cell) {
	if (capacity>0 && size() == capacity) throw CellBufferOverflow();
	Clist.push_front(cell);
}

inline Cell* CellBintree::pop() {
	Cell* c = Clist.front();
	Clist.pop_front();
	return c;
}

inline Cell* CellBintree::top() const {
	if(Clist.empty())
		return NULL;
	else
		return Clist.front();
};


} // end namespace ibex
#endif // __IBEX_CELL_BINTREE_H__

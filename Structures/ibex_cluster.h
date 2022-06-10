

//============================================================================
//                                  I B E X
// File        : ibexCellBdfs.h
// Author      : Thomas Richard de Latour
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : June, 00
// Last Update : June, 00
//============================================================================
#ifndef __IBEX_CLUSTER_H__
#define __IBEX_CLUSTER_H__

#include "ibex.h"
#include "Logger.h" // Logger to store and print data in csv

namespace ibex {

struct CellCluster{

	IntervalVector Iref;

	std::list<IntervalVector> clist;

	std::vector<double> vdist;

	int id;

};


/** \ingroup strategy
 *
 * \brief Cluster grouping cells for clustering algorithm
 */
class Cluster {
public:


	Cluster(std::vector<CellCluster> cl, BitSet bit, IntervalVector init);

	Cluster(std::list<Cell*> cl, BitSet bit, IntervalVector init);

	Cluster(std::list<IntervalVector> boxes, BitSet bit, IntervalVector init);


	void init_matrix();

	void group_min();

	void update_matrix();

	void find_min();

	int size();

	void print();

	void write_all(System sys);

protected :

	pair<int,int> min_clust;

	std::vector<CellCluster> clustlist;

	double distance_clust(IntervalVector x1, IntervalVector x2);

	BitSet bitset;

	IntervalVector init_box;

};


inline Cluster::Cluster(std::vector<CellCluster> cl, BitSet bit, IntervalVector init): clustlist(cl), bitset(bit), init_box(init){

}

inline Cluster::Cluster(std::list<Cell*> cl, BitSet bit, IntervalVector init):  bitset(bit), init_box(init){
	int k = 0;
	for(Cell* cell : cl){
		CellCluster clust;
		clust.Iref = cell->box;
		clust.clist.push_back(cell->box);
		clust.id = k;
		clust.vdist.resize(cl.size());
		clustlist.push_back(clust);
		k+=1;
	}

}

inline Cluster::Cluster(std::list<IntervalVector> boxes, BitSet bit, IntervalVector init): bitset(bit), init_box(init){
	int k = 0;
	for(IntervalVector box : boxes){
		CellCluster clust;
		clust.Iref = box;
		clust.clist.push_back(box);
		clust.id = k;
		clust.vdist.resize(boxes.size());
		clustlist.push_back(clust);
		k+=1;
	}

}

inline int Cluster::size(){
	return clustlist.size();
}

inline void Cluster::init_matrix(){ // initialise les distances entre les clusters, le regroupement des indices
	double mind = POS_INFINITY;
	double dist;
	min_clust.first = 0;
	min_clust.second = 0;
	for(int i =0; i<clustlist.size() ; i++){
		for(int j = i+1; j<clustlist.size() ; j++){
			dist = distance_clust(clustlist[i].Iref[bitset],clustlist[j].Iref[bitset]);
			clustlist[i].vdist[j] = dist;
			if(dist<mind){
				mind =dist;
				min_clust.first = i;
				min_clust.second = j;
			}
		}
	}
}

inline void Cluster::group_min(){
	clustlist[min_clust.first].Iref |=  clustlist[min_clust.second].Iref;
	for(IntervalVector box : clustlist[min_clust.second].clist){ // regroupe les deux cluster les plus proches à l'indice min_clust.first
		clustlist[min_clust.first].clist.push_back(box);
	}	// min_clust.first = indice du dernier regroupement
	clustlist.erase(clustlist.begin()+min_clust.second);	// elimine le second cluster
	for(CellCluster cc : clustlist){	// adapte la distance des cluster
		cc.vdist.erase(cc.vdist.begin()+min_clust.second);
		if(cc.id > min_clust.second) cc.id = cc.id-1 ;
	}	// min_clust.second = indice du cluster éliminer pour le regroupement
}

inline void Cluster::update_matrix(){ // min_clust.first = indice du dernier regroupement
	for(int i = 0; i < min_clust.first; i++){
		clustlist[i].vdist[min_clust.first]=distance_clust(clustlist[i].Iref[bitset],clustlist[min_clust.first].Iref[bitset]);
	}
	for(int i =min_clust.first+1 ; i< clustlist.size();i++){
		clustlist[min_clust.first].vdist[i] = distance_clust(clustlist[i].Iref[bitset],clustlist[min_clust.first].Iref[bitset]);
	}
}

inline void Cluster::find_min(){
	double mind = POS_INFINITY;
	double dist;
	for(int i = 0; i < clustlist.size() ; i++){
		for(int j= i+1; j < clustlist.size() ; j++){
			dist = clustlist[i].vdist[j];
			if(mind > dist){
				mind =dist;
				min_clust.first = i;
				min_clust.second = j;
			}
		}
	}
}

inline double Cluster::distance_clust(IntervalVector x1, IntervalVector x2){
	assert((x1.size()==x2.size()) && (x1.size()==init_box.size()));

	IntervalVector crit_box = init_box[bitset];

	double max = distance(x1[0],x2[0])/crit_box[0].diam();
	for (int i=1; i<x1.size(); i++) {
		if (max<distance(x1[i],x2[i])/crit_box[i].diam()) max = distance(x1[i],x2[i])/crit_box[i].diam();
	}
	return max;
}


inline void Cluster::print(){
	cout <<"CLUSTER"<<endl;
	for(auto c : clustlist){
		cout <<"\t" << c.id ;
		cout << " "<<c.clist.size();
		cout << " "<<c.Iref;
		cout<<endl;
	}
	cout<<"/n"<<endl;

}

inline void Cluster::write_all(System sys){
	string filename;
	cin >> filename;
	ofstream ofs (filename, ofstream::out);

	for(int p =0 ; p <sys.args.size(); p++){
		ofs <<";;;"<< sys.args[p].name;
	}
	ofs <<endl;

	for(int i =0; i<clustlist.size(); i++){
		ofs <<"Cluster Number : ;"<<i<<";center : ;";
		for(int j=0; j<clustlist[i].Iref.size();j++){
			ofs << clustlist[i].Iref[j].lb() << ";" << clustlist[i].Iref[j].ub();
			ofs << ";;";
		}
		ofs <<endl;
		for(auto box : clustlist[i].clist){
			ofs<<";";
			for(int k=0; k<box.size();k++){
						ofs <<";;"<< box[k].lb() << ";" <<box[k].ub();
					}
			ofs <<endl;

		}
	}
	ofs.close();
}
} // end namespace ibex
#endif // __IBEX_CELL_CLUSTER_H__


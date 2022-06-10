//============================================================================
//                                  I B E X
// File        : Integer First bisector
// Author      : Thomas Richard de Latour
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : July 2020
//============================================================================
#ifndef LOGGER_H_
#define LOGGER_H_

#include <math.h>
#include "ibex.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>

using namespace std;

// Logger to save exploration information when using Anytime BPA Algorithm


namespace ibex {

//Normalize hausdorff distance between two boxes on their initial domain
double distance_norm(const IntervalVector& x1, const IntervalVector& x2,const IntervalVector& init_box) {
	assert((x1.size()==x2.size()) && (x1.size()==init_box.size()));

	double max = distance(x1[0],x2[0])/init_box[0].diam();
	for (int i=1; i<x1.size(); i++) {
		if (max<distance(x1[i],x2[i])/init_box[i].diam()) max = distance(x1[i],x2[i])/init_box[i].diam();
	}
	return max;
}

// Volume of an interval vector
double Intvolume(const IntervalVector& x1, Vector prec) {
	assert((x1.size()==prec.size()));
	double vol = 1.0;
	for(int i=0; i<x1.size(); i++){
		if((!x1[i].is_degenerated()) && (x1[i].diam() >= prec[i])){
			vol = vol*x1[i].diam()/prec[i];
		}
	}
	return vol;
}



class Logger : public CovList {
public:
	/**
	 * \brief Create a new, empty covering list.
	 *
	 * \param n - the dimension of the covered set.
	 */
	Logger(size_t n);

	/**
	 * \brief Conversion from a CovList.
	 */
	Logger(const CovList& cov, const vector<double> covTime,const vector<int> covIt,const vector<int> covUnex,const vector<int> covSplit);

	/**
	 * \brief Delete this
	 */
	~Logger();
	/**
	 * \brief Number of boxes
	 */
	size_t size() const;

	/**
	 * \brief add a new box
	 */

	void add(const IntervalVector& x);
	/**
	 * \brief Get the ith box.
	 */
	const IntervalVector& operator[](int i) const;

	/*	 * brief write only the main information	 */
	void write_csv(const string filename);

	void write_all(const string path,int nb);

	void write_txt(const string path,int nb);

	void write_resume(const string file);

	/* * brief name of the extracting strategy */
	string Strat;

	/* * brief name of problem */
	string Sys;

	/**
	 * \brief Conversion from a CovList.
	 */
	void copy(const CovList& cov);

	/**
	 * \brief Conversion from data
	 */
	void set_data(const vector<double> covTime,const vector<int> covIt,const vector<int> covUnex,const vector<int> covSplit);

	/**
	 * \brief reset the CovList except the hypermesh*/
	void reset();
	
	void set_box(IntervalVector& box);

	void set_box(IntervalVector& box, std::vector<std::string> names, BitSet bitset);
	
	bool norm;

	Vector* prec;

	std::pair<double,double> Qfinal;

	double get_QA();

	double get_QH();

	std::list<IntervalVector> get_solutions(); // return the list of the nb first solution
	
	IntervalVector* init_box ;
	
	BitSet var;
	
protected:

	/* * brief Average minimum distance Quality indicator	 */
	void QA();


	/* * brief enveloppe Quality indicator	 */
	void QH();


	struct Data {
		std::list<IntervalVector> lst;
		std::vector<IntervalVector*> vec;
		std::vector<double> time;
		std::vector<int> it;
		std::vector<int> unex;
		std::vector<int> split;
	} *data;

	struct Qindic {
		std::vector<double> QA;
		std::vector<double> QH;
	} *qindic;

	bool own_data;
	
	std::vector<std::string> var_names;
	
};
/*================================== inline implementations ========================================*/
inline Logger::Logger(size_t n): CovList(n), data(new Data()), qindic(new Qindic()),own_data(true){
	norm = false;
	init_box = NULL;
}

inline Logger::Logger(const CovList& cov, const vector<double> covTime,const vector<int> covIt,const vector<int> covUnex,const vector<int> covSplit):CovList(cov.size()){
	const Logger* covlist = dynamic_cast<const Logger*>(&cov);
	norm = false;
	if (covlist) {
		data = new Data();
		qindic =new Qindic();
		data->lst  = covlist->data->lst;
		for (list<IntervalVector>::iterator it=data->lst.begin(); it!=data->lst.end(); ++it)
			data->vec.push_back(&(*it));
		own_data = true;
		data->time = covTime;
		data->it = covIt;
		data->unex = covUnex;
		data->split = covSplit;
		init_box = NULL;
		prec = NULL;

	} else {
		data = new Data();
		qindic =new Qindic();
		// nothing to add
		own_data = true;
	}
}

inline void Logger::copy(const CovList& cov){
	this->reset();

	for (int i=0;i<cov.size();i++)
		(*this).add(cov[i]);

}

inline void Logger::set_data(const vector<double> covTime,const vector<int> covIt,const vector<int> covUnex,const vector<int> covSplit){
	data->time = covTime;
	data->it = covIt;
	data->unex = covUnex;
	data->split = covSplit;
}

inline void Logger::reset(){
	data->lst.clear();
	data->vec.clear();
	data->time.clear();
	data->it.clear();
	data->unex.clear();
	data->split.clear();
	qindic->QA.clear();
	qindic->QH.clear();
	//init_box = NULL;
	//prec = NULL;
}

inline size_t Logger::size() const {
	return data->lst.size();
}

inline void Logger::add(const IntervalVector& x) {

	if (n!=(size_t) x.size())
		ibex_error("[CovList] boxes must have all the same size.");
	data->lst.push_back(x);
	data->vec.push_back(&data->lst.back());
}


inline const IntervalVector& Logger::operator[](int i) const {
	return *(data->vec)[i];
}

inline Logger::~Logger() {
	if (own_data) {
		delete data;
		delete qindic;
	}
}

//average min hausdorff distance quality indicator computation
inline void Logger::QA(){
	qindic->QA.push_back(0);
	if(this->size()>1){
		for (int j =1;j<this->size(); j++){
			double Q = qindic->QA.back()*j;
			double d = POS_INFINITY;
			if(norm){
				for (int i =0;i<j; i++){
					if(distance_norm((*this)[j][var],(*this)[i][var],(*init_box)[var])<d){
						d = distance_norm((*this)[j][var],(*this)[i][var],(*init_box)[var]);
					}
				}
			}else{
				for (int i =0;i<j; i++){
					if(distance((*this)[j][var],(*this)[i][var])<d){
						d = distance((*this)[j][var],(*this)[i][var]);
					}
				}
			}
			Q += d;
			qindic->QA.push_back(Q/(j+1));
		}
	}
}

// Hypervolume quality indicator computation
inline void Logger::QH(){
	IntervalVector box;
	if(this->size()>0){
		box = (*this)[0];
		qindic->QH.push_back(Intvolume(box[var],*prec)/Intvolume((*init_box)[var],*prec));
		for(int i=1; i<this->size();i++){
			box = box|(*this)[i];
			qindic->QH.push_back(Intvolume(box[var],*prec)/Intvolume((*init_box)[var],*prec));
		}
	}else{
		box = IntervalVector(1);
		qindic->QH.push_back(box.volume()/Intvolume((*init_box)[var],*prec));
	}
}

// gete quality indicators
inline double Logger::get_QA(){
	if(qindic->QA.empty())
		this->QA();
	return qindic->QA[this->size()-1];
};

inline double Logger::get_QH(){
	if(qindic->QH.empty())
		this->QH();
	return qindic->QH[this->size()-1];
};

inline void Logger::set_box(IntervalVector& box){
	init_box = &box;
	var = var.all(init_box->size());
	for(int i=0; i< init_box->size();i++){
	}
}

inline void Logger::set_box(IntervalVector& box,std::vector<std::string> names, BitSet bitset){
	init_box = &box;
	var = bitset;
	var_names = names;
	for(int i=0; i< init_box->size();i++){
	}
}

// Save computation info for each solution in a row of a csv file
inline void Logger::write_csv(const string path){
	string filename=path+"QI_"+Sys+".csv";
	ofstream ofs;
	ofs.open(filename,ios::out | ios::app);

	this->QA();
	this->QH();	
	for(int i =0; i<this->size(); i++){
		ofs << Strat << ",";
		ofs << data->time[i] <<",";
		ofs << i+1 <<",";
		ofs << qindic->QA[i] <<",";
		ofs << qindic->QH[i] <<",";
		ofs << data->split[i] <<",";
		ofs << data->it[i] <<",";
		ofs << data->unex[i];
		ofs <<endl;
	}

	ofs.close();

}

// Save final computation info in a row of a csv file
inline void Logger::write_resume(const string path){
	string filename=path+"RESUME_"+Sys+".csv";
	ofstream ofs;
	ofs.open(filename,ios::out | ios::app);

	this->QA();
	this->QH();
	int i =this->size()-1;
	ofs << Sys << ",";
	ofs << Strat << ",";
	ofs << data->time[i] <<",";
	ofs << i+1 <<",";
	ofs << qindic->QA[i] <<",";
	ofs << qindic->QH[i] <<",";
	ofs << data->split[i] <<",";
	ofs << data->it[i] <<",";
	ofs << data->unex[i];
	ofs <<endl;
	ofs.close();

}

// CSV with all computation details
inline void Logger::write_all(const string path, int nb){
	string filename=path+"ALL_"+Sys+"_"+Strat+".csv";//round(CPU_Time*1000)
	ofstream ofs (filename, ofstream::out);
	ofs <<endl;
	ofs << "Strategy;" << Strat<<endl;
	ofs << "Nb_solutions;" << this->size()<<endl;
	ofs << "Unexplored_boxes;" << data->unex[-1]<<endl;
	ofs << "CPU_Time;" << data->time[nb]<<endl;
	ofs << "Nb_Splits;" << data->split[nb]<<endl;
	ofs << "Nb_iterations;" << data->it[nb]<<endl;

	for(int i =0; i<nb; i++){
		for(int j=0; j<data->vec[0]->size();j++){
			ofs << data->vec[i][0][j].lb() << ";" << data->vec[i][0][j].ub();
			ofs << "  ";
		}
		ofs <<endl;
	}
	ofs.close();
}

// CSV with all solutions saved as middle point of the solution boxes
inline void Logger::write_txt(const string path,int nb){
	string filename=path+"MDP_"+Sys+"_"+Strat+".csv";//round(CPU_Time*1000)
	ofstream ofs (filename, ofstream::out);

	for(int j=0; j<data->vec[0]->size();j++){
		ofs << var_names[j] << ",";
	}
	ofs <<endl;
	for(int i =0; i<nb; i++){
		for(int j=0; j<data->vec[0]->size();j++){
			ofs << (data->vec[i][0][j].lb()+data->vec[i][0][j].ub())/2 << ",";
		}
		ofs <<endl;
	}
	ofs.close();
}

inline std::list<IntervalVector> Logger::get_solutions(){
	return data->lst;
}

}

#endif /* LOGGER_H_ */

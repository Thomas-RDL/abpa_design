#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <typeinfo> 
#include <string>
#include <filesystem>
#include <dirent.h>
#include <algorithm>

#include "Structures/Anytime_BPA.h" // Anytime Branch and Prune Algorithm + criteria function strategies
#include "Structures/Logger.h" // Logger to store and print data in csv
#include "Structures/ibex_LFNorm.h" // Normalized Largest First
#include "Structures/ibex.h"
#include "Structures/ibex_table_ctc.h"

using namespace std;
using namespace ibex;
#define stringer( x ) printf(#x)

/*
 * main for Anytime BPA_BFS
 *  Created on: 28 oct. 2020
 *      Author: Thomas Richard de Latour
 */
 */

 
int main(){

	std::vector<std::pair<std::string, std::vector<double>>> cell = read_csv("cells_data.csv");

	// var for ctc_table
	std::vector<std::string> vars ={"n","Cc","Uc","Tc","Wc","Hc","Mc","Type","Form","Cprice","Rint","Ucmin","Crate_peak"};
	int k =0;

	Matrix data_cell(vars.size(),cell.at(0).second.size());

	for(int i=0;i < cell.size();i++){
		if(std::find(vars.begin(), vars.end(), cell.at(i).first) != vars.end()) {
			for(int j=0;j < data_cell.nb_cols();j++){
				data_cell[k][j] = cell.at(i).second.at(j);
			}
			k+=1;
		}
	}
	
	cout << "***** CELL DATA ****"<< endl;

	cout<<data_cell<<endl;

	string path ="./Results/"; //path for results
	string sysName ="battery"; //problem sysName

	cout << "*************** Problem name : "<< sysName<<" **************"<< endl;

//	string syscar = "/home/E19H844K/Documents/these-thomas-ls2n/Articles/EAAI/Code/models/"+sysName+".mbx";
	string syscar = sysName+".mbx";
	const char *syschar = syscar.c_str();
	System sys(syschar); // problem path
	cout << sys <<endl;
	
	sysName ="battery_gwp_rev_cost"; // result files prefix name
	
	//Resolution details
	//  VAR :   n, Cc, Uc, Tc, Wc, Hc, Mc, Type, Form, Cprice, Rint, Ucmin, Crate_peak, Ns, Np, Ubp, Ebp, Nsy, M_batt, V_batt, Ppeak, GWP, REV, Cost, E_batt, E_ri
	// Design variables 	
	BitSet bset(sys.nb_var);	
	list<int> indice_var{0,13,14,21,22,23}; //  correspond to variables n, Ns, Np, GWP, REV, Cost
	for(auto i = indice_var.begin();i!=indice_var.end();i++){
		bset.add(*i);
	}

	// Precision of the computation for each variable
	double eps=1e-3; // mean precision
	std::initializer_list<double> l_prec{1,0.1,0.01,0.00001,0.00001,0.00001,0.001,1,1,0.1,0.01,0.01,0.1,1,1,1,1,1,0.1,0.001,1,1,0.1,1,1,1};
	Vector V_eps(l_prec);// vector of precision detailed
	assert(l_prec.size() = sys.nb_var);	
	int max_sol = 2500; // Maximum number of solutions
	IntervalVector init_box(sys.box); // initial box of the system

	//*************Bisector (utile = bsc)**************
	LargestFirst LF(V_eps,0.5); // Bissector Largest First
	RoundRobin RR(V_eps,0.5); // Bissector Round Robin
	LFNorm LFN(V_eps,0.5); // Bissector Largest First Normalize with initial domain and precision
	BscFirst BFisrt(bset,RR, RR,V_eps,0.5); // Bisector which priorize bissection on selected dimenssions
	SmearSumRelative SF(sys,V_eps); // Smear Sum Bissector

	Bsc* bsc = &RR;

	//******** Contractor (utile = ctc)************
	CtcTable ctc_tab(sys,data_cell,vars);

	// HC4
	CtcHC4 ctc_hc4(sys); 
	CtcFixPoint ctc_fphc4(ctc_hc4);

	// 3BCID
	Ctc3BCid ctc_3bcid(ctc_hc4);

	// ACID
	CtcAcid ctc_acid(sys, ctc_hc4);
	
	// Contractor composition
	CtcCompo ctcvar(ctc_fphc4,ctc_acid,ctc_tab);
	
	// Contractor fixpoint
	CtcFixPoint ctc_bis(ctcvar); 
	ctc_bis.contract(init_box); // first contraction to avoid infinite bounds (application case)
	
	// Print first contraction
	cout << init_box<<endl;
	 

	//****Buffers******* 
	CellBintree buff_bt;// BestFS multiset + stack
	CellStack buff_dfs; // Depth First Search
	CellList buff_breadth; // Breadth First Search
	CellWidth buff_width; // Largest First search 

	// Set the initial box of the problem
	LFN.set_box(init_box);  

	Logger* S =new Logger(sys.nb_var);
	init_CSV(path,"RESUME_"+sysName,eps,typeid(RR).name(),typeid(ctc_bis).name(),false); // headers for resume result
	init_CSV(path,"QI_"+sysName,eps,typeid(RR).name(),typeid(ctc_bis).name(),false); // headers for quality indicators result
	
	//Init Logger  
	S->norm = true; // normalization fo QA and QH 
	S->Sys = sysName;
	S->set_box(init_box, sys.var_names(), bset);
	S->prec = &V_eps;

	//Depth-First exploration
	//anytime_BFS_DFS(ctc_bis,*bsc,buff_bt,*S,crit_idfs,30,10000,false);
	
	// IDFS exploration strategy
	S->Strat="_idfs_"+to_string(max_sol);
	anytime_BFS_DFS(ctc_bis,*bsc,buff_bt,*S,crit_idfs,timeout,max_sol,true);
	S->write_resume(path);
	S->write_txt(path, S->size()); 
	S->write_csv(path);
	

	// Normalized DMDFS exploration strategy
	buff_bt.dynamic =true;	
	S->Strat="_dmdfs_norm_"+to_string(max_sol);
	anytime_BFS_DFS(ctc_bis,*bsc,buff_bt,*S,crit_dmdfs_norm,timeout,max_sol,true);
	S->write_resume(path);
	S->write_txt(path, S->size());
	S->write_csv(path);
	
	buff_bt.dynamic =false; 
	

	cout <<"\t ----> FIN"<<endl;
	return 0; 
}



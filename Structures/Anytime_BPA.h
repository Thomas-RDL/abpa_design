/*
 * Anytime_BPA.h
 *
 *  Created on: 28 oct. 2020
 *      Author: Thomas Richard de Latour
 */

#include <fstream>
#include "ibex_CellBintree.h"
#include <iostream>

#include "ibex.h"
#include "Logger.h"
using namespace std;
using namespace ibex;
#define stringer( x ) printf(#x)

bool is_too_small(const IntervalVector& box, double eps){
	for (int i=0; i<box.size(); i++)
		if (box[i].diam()>eps) return false;
	return true;
}

bool V_is_too_small(const IntervalVector& box, Vector eps){
	assert(box.size()==eps.size());
	for (int i=0; i<box.size(); i++)
		if (box[i].diam()>eps[i]) return false;
	return true;
}

/*
 * Measure-of-best functions
 */

/*/
 * No criteria
 */

double no_crit(const Logger& cov, CellBuffer& buffer, Cell& cell){return 0.0;}

/*/
 * Measure of best function for IDFS
 */

double crit_idfs(const Logger& cov, CellBuffer& buffer, Cell& cell){
	double crit = cell.depth;
	//cout << cell.depth<<endl;
	return -crit;
}

/*/
 * Measure of best function for RDFS
 */
double crit_rand(const Logger& cov, CellBuffer& buffer, Cell& cell){
    double crit=1.0*std::rand()/RAND_MAX;
    return crit;
}

/*/
 * Measure of best function for DLFS
 */

double crit_ldfs(const Logger& cov, CellBuffer& buffer, Cell& cell){
	return cell.box.max_diam();
}

double crit_ldfs_norm(const Logger& cov, CellBuffer& buffer, Cell& cell){
	double crit = 0;
	for(int i=0; i< cell.box.size();i++ ){
		if(cell.box[i].diam()/((*cov.init_box)[i].diam()) > crit)
				crit = cell.box[i].diam()/((*cov.init_box)[i].diam());
	}
	return crit;
}

double crit_width_bitset(const Logger& cov, CellBuffer& buffer, Cell& cell){
	double crit = 0;
	for(int i=0; i< cell.box.size();i++ ){
		if(cell.box[i].diam()/((*cov.init_box)[i].diam()) > crit && cov.var[i])
				crit = cell.box[i].diam()/((*cov.init_box)[i].diam());
	}
	return crit;
}

double crit_perimeter(const Logger& cov, CellBuffer& buffer, Cell& cell){
	return cell.box.perimeter();
}

double crit_vol(const Logger& cov, CellBuffer& buffer, Cell& cell){
	return cell.box.volume();
}

/*/
 * Measure of best function for DMDFS
 */

double crit_dmdfs(const Logger& cov, CellBuffer& buffer, Cell& cell){
	//critère de tri
	double crit = POS_INFINITY;
	if(cov.size()>=1){
		crit = distance(cell.box,cov[0]);
		for (int i=1 ; i < cov.size() ; i++){
			double d = distance(cell.box,cov[i]);
			if(d<crit) crit = d;
		}
	}
	return crit;
}

double crit_dmdfs_norm(const Logger& cov, CellBuffer& buffer, Cell& cell){
	//critère de tri
	double crit = POS_INFINITY;
	if(cov.size()>=1){
			crit = distance_norm(cell.box,cov[0],*cov.init_box);
			for (int i=1 ; i < cov.size() ; i++){
				double d = distance_norm(cell.box,cov[i],*cov.init_box);
				if(d<crit) crit = d;
			}

	}

	return crit;
}

double crit_dmdfs_bitset(const Logger& cov, CellBuffer& buffer, Cell& cell){
	//critère de tri
	double crit = POS_INFINITY;
	if(cov.size()>=1){
			crit = distance_norm(cell.box[cov.var],cov[0][cov.var],(*cov.init_box)[cov.var]);
			for (int i=1 ; i < cov.size() ; i++){
				double d = distance_norm(cell.box[cov.var],cov[i][cov.var],(*cov.init_box)[cov.var]);
				if(d<crit) crit = d;
			}

	}

	return crit;
}

/*
 * Initiate a .csv file which will contain informations on the resolution before post-evaluation
 */

void init_CSV(const string path, const string file,double prec,string bsc,string ctc, bool zero){
	string filename=path+file+".csv";
	ofstream ofs;
	if(zero){
		ofs.open(filename,ios::out);
	}else{
		ofs.open(filename,ios::out | ios::app);
	}
	if(ofs.tellp()==0){
//		ofs << "file_name,precision,bissecteur,contracteur" <<endl;
//		ofs << file <<","<<prec<<","<<bsc<<","<<ctc<<endl;
		ofs << "Pbm,CPU_Time,Number_Solutions,QA,QH,Number_Split,Number_It,Unexplored_Boxes" <<endl;
		ofs.close();
	}
}

/*
 * Anytime branch and prune BFS with data logger
 */

void anytime_BFS_DFS( Ctc& ctc, Bsc& bsc, CellBuffer& L,Logger& Sat, double (*critere)(const Logger&, CellBuffer&, Cell&), double timeout, int max_sol, bool QI){

	Timer time_limit;
	time_limit.start();

	L.flush(); // unexplored boxes list

	Cell* root=new Cell(*Sat.init_box); // initial box
	L.push(root);

	Sat.reset();

	int nb_it = 0; // number of iteration
	int nb_split= 0; //number of bissection
	int nb_fail = 0;//
	int nb_Sol = 0; // number of epsilon-boxes
	vector<double> covTime;
	vector<int> covSplit;
	vector<int> covIt;
	vector<int> covUnex;

	cout<<">>>> start - ";

	while(!L.empty() && nb_Sol < max_sol){ // Write the intermediates results in a global .csv

		CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
		if (sb!=NULL){
			sb->init_list();
		}
		
		//cout << " "<<L.size()<<endl;

		try {
			time_limit.check(timeout);
		}
		catch(TimeOutException& e) {
			break;
		}

		Cell* cell= L.pop();
		nb_it++;
		ctc.contract(cell->box);
		if(!cell->box.is_empty()){
			if(!V_is_too_small(cell->box,*Sat.prec)){

				pair<Cell*,Cell*> new_cells=bsc.bisect(*cell);
				nb_split++;

				L.push(new_cells.first);
				L.push(new_cells.second);



			}
			else {

				Sat.add(cell->box);
				covTime.push_back(time_limit.get_time());
				covIt.push_back(nb_it);
				covSplit.push_back(nb_split);
				covUnex.push_back(L.size());

				nb_Sol++;

				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
				if (sb!=NULL){
					sb->update(Sat,(*critere));
				}

			}
		}else{
			nb_fail++;
			if(nb_fail>1000){
				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
					if (sb!=NULL){
						sb->update(Sat,(*critere));
					}
					
			};
		
		}
		delete cell;
	}

	/*
	 * store data results in a Logger object (see Logger.h)
	 */

	time_limit.stop();

	Sat.set_data(covTime,covIt,covUnex,covSplit);

	// print results
	cout << "-- Stratégie : " << Sat.Strat << endl;
	cout << "\t Temps de résolution : " << time_limit.get_time() << endl;
	cout << "\t Nombre de solutions : " << Sat.size()<<endl<<endl;
	
	CellBintree* sb = dynamic_cast<CellBintree*>(&L);
	if (sb!=NULL)cout << "\t Pending list       : " <<sb->size_tree()<<endl<<endl;
	
	cout << "\t nombre de split       : " << nb_split<<endl<<endl;
	cout << "\t nombre de itération       : " << nb_it<<endl<<endl;
	
	if(QI){
		cout << "\t Indicateur QA       : " << Sat.get_QA()<<endl<<endl;
		cout << "\t Indicateur QH       : " << Sat.get_QH()<<endl<<endl;
		}

}

void anytime_BFS_DFS_2( Ctc& ctc, Bsc& bsc, CellBuffer& L,Logger& Sat, double (*critere)(const Logger&, CellBuffer&, Cell&), double timeout, int max_sol){

	Timer time_limit;
	time_limit.start();

	L.flush(); // unexplored boxes list

	Cell* root=new Cell(*Sat.init_box); // initial box
	L.push(root);

	Sat.reset();

	int nb_it = 0; // number of iteration
	int nb_split= 0; //number of bissection
	int nb_fail = 0;//
	int nb_Sol = 0; // number of epsilon-boxes
	vector<double> covTime;
	vector<int> covSplit;
	vector<int> covIt;
	vector<int> covUnex;

	cout<<">>>> start - ";

	while(!L.empty() && nb_Sol < max_sol){ // Write the intermediates results in a global .csv

		CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
		if (sb!=NULL){
			sb->init_list();
		}
		
		//cout << " "<<L.size()<<endl;

		try {
			time_limit.check(timeout);
		}
		catch(TimeOutException& e) {
			break;
		}

		Cell* cell= L.pop();
		nb_it++;
		ctc.contract(cell->box);
		if(!cell->box.is_empty()){
			if(!V_is_too_small(cell->box,*Sat.prec)){

				pair<Cell*,Cell*> new_cells=bsc.bisect(*cell);
				nb_split++;

				L.push(new_cells.first);
				L.push(new_cells.second);



			}
			else {

				Sat.add(cell->box);
				covTime.push_back(time_limit.get_time());
				covIt.push_back(nb_it);
				covSplit.push_back(nb_split);
				covUnex.push_back(L.size());

				nb_Sol++;

				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
				if (sb!=NULL){
					sb->update(Sat,(*critere));
				}

			}
		}
		delete cell;
	}

	/*
	 * store data results in a Logger object (see Logger.h)
	 */

	time_limit.stop();

	Sat.set_data(covTime,covIt,covUnex,covSplit);

	// print results
	cout << "-- Stratégie : " << Sat.Strat << endl;
	cout << "\t Temps de résolution : " << time_limit.get_time() << endl;
	cout << "\t Nombre de solutions : " << Sat.size()<<endl<<endl;
	
	CellBintree* sb = dynamic_cast<CellBintree*>(&L);
	if (sb!=NULL)cout << "\t Pending list       : " <<sb->size_tree()<<endl<<endl;
	
	cout << "\t nombre de split       : " << nb_split<<endl<<endl;
	cout << "\t nombre de itération       : " << nb_it<<endl<<endl;
	cout << "\t Indicateur QA       : " << Sat.get_QA()<<endl<<endl;
	cout << "\t Indicateur QH       : " << Sat.get_QH()<<endl<<endl;

}

void anytime_MDFS( Ctc& ctc, Bsc& bsc, CellBuffer& L,Logger& Sat, double (*critere)(const Logger&, CellBuffer&, Cell&), double timeout, int max_sol){

	Timer time_limit;
	time_limit.start();

	L.flush(); // unexplored boxes list

	Cell* root=new Cell(*Sat.init_box); // initial box
	L.push(root);

	Sat.reset();

	int nb_it = 0; // number of iteration
	int nb_split= 0; //number of bissection
	int nb_Sol = 0; // number of epsilon-boxes
	vector<double> covTime;
	vector<int> covSplit;
	vector<int> covIt;
	vector<int> covUnex;

	cout<<">>>> start - ";

	while(!L.empty() && nb_Sol < max_sol){ // Write the intermediates results in a global .csv

		CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
		if (sb!=NULL){
			sb->init_list();
		}


		try {
			time_limit.check(timeout);
		}
		catch(TimeOutException& e) {
			break;
		}
		Cell* cell= L.pop();
		nb_it++;
		ctc.contract(cell->box);
		if(!cell->box.is_empty()){
			if(!V_is_too_small(cell->box,*Sat.prec)){

				pair<Cell*,Cell*> new_cells=bsc.bisect(*cell);
				nb_split++;

				L.push(new_cells.first);
				L.push(new_cells.second);

				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
				if (sb!=NULL){
					sb->insert(Sat,(*critere));
				}

			}
			else {
				Sat.add(cell->box);

				covTime.push_back(time_limit.get_time());
				covIt.push_back(nb_it);
				covSplit.push_back(nb_split);
				covUnex.push_back(L.size());

				nb_Sol++;

				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
				if (sb!=NULL){
					sb->update(Sat,(*critere));
				}

			}
		}
		delete cell;
	}

	/*
	 * store data results in a Logger object (see Logger.h)
	 */

	time_limit.stop();


	Sat.set_data(covTime,covIt,covUnex,covSplit);

	// print results
	cout << "-- Stratégie : " << Sat.Strat << endl;
	cout << "\t Temps de résolution : " << time_limit.get_time() << endl;
	cout << "\t Nombre de solutions : " << Sat.size()<<endl<<endl;
//	cout << "\t Indicateur QA       : " << Sat.get_QA()<<endl<<endl;
//	cout << "\t Indicateur QH       : " << Sat.get_QH()<<endl<<endl;

}

std::vector<std::pair<std::string, std::vector<double>>> read_csv(std::string filename){
	// Reads a CSV file into a vector of <string, vector<int>> pairs where
	// each pair represents <column name, column values>

	// Create a vector of <string, int vector> pairs to store the result
	std::vector<std::pair<std::string, std::vector<double>>> result;

	// Create an input filestream
	std::ifstream myFile(filename);

	// Make sure the file is open
	if(!myFile.is_open()) throw std::runtime_error("Could not open file");

	// Helper vars
	std::string line, colname;
	double val;

	// Read the column names
	if(myFile.good())
	{
		// Extract the first line in the file
		std::getline(myFile, line);

		// Create a stringstream from line
		std::stringstream ss(line);

		// Extract each column name
		while(std::getline(ss, colname, ',')){

			// Initialize and add <colname, int vector> pairs to result
			result.push_back({colname, std::vector<double> {}});
		}
	}

	// Read data, line by line
	while(std::getline(myFile, line))
	{
		// Create a stringstream of the current line
		std::stringstream ss(line);

		// Keep track of the current column index
		int colIdx = 0;

		// Extract each integer
		while(ss >> val){

			// Add the current integer to the 'colIdx' column's values vector
			result.at(colIdx).second.push_back(val);

			// If the next token is a comma, ignore it and move on
			if(ss.peek() == ',') ss.ignore();

			// Increment the column index
			colIdx++;
		}
	}

	// Close file
	myFile.close();

	return result;
}

//void anytime_BFS_DFS(const string pathname, Ctc& ctc, Bsc& bsc, CellBuffer& L,Logger& Sat, IntervalVector& init_box,double (*critere)(const CovList&, CellBuffer&, Cell&), Vector epsilon, double timeout, int max_sol, bool dynamic){
//
//	Timer time_limit;
//	time_limit.start();
//
//	L.flush(); // unexplored boxes list
//
//	IntervalVector box(init_box);
//
//	Cell* root=new Cell(box); // initial box
//	L.push(root);
//
//	CovList S(init_box.size());
//	Sat.reset();
//	Sat.set_box(init_box);
//
//	int nb_it = 0; // number of iteration
//	int nb_split= 0; //number of bissection
//	int nb_Sol = 0; // number of epsilon-boxes
//	vector<double> covTime;
//	vector<int> covSplit;
//	vector<int> covIt;
//	vector<int> covUnex;
//
//	cout<<">>>> start - ";
//
//	while(!L.empty() && nb_Sol < max_sol){ // Write the intermediates results in a global .csv
//
//		CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
//		if (sb!=NULL){
//			sb->init_list();
//		}
//
//
//		try {
//			time_limit.check(timeout);
//		}
//		catch(TimeOutException& e) {
//			break;
//		}
//
//		Cell* cell= L.pop();
//
//
//		nb_it++;
//		ctc.contract(cell->box);
//		if(!cell->box.is_empty()){
//			if(!V_is_too_small(cell->box,epsilon)){
//
//				pair<Cell*,Cell*> new_cells=bsc.bisect(*cell);
//				nb_split++;
//
//				L.push(new_cells.first);
//				L.push(new_cells.second);
//
//			}
//			else {
//
//				S.add(cell->box);
//
//				covTime.push_back(time_limit.get_time());
//				covIt.push_back(nb_it);
//				covSplit.push_back(nb_split);
//				covUnex.push_back(L.size());
//
//				nb_Sol++;
//
//				CellBintree* sb = dynamic_cast<CellBintree*>(&L); // case multiset + list
//				if (sb!=NULL){
//					sb->update(S,(*critere));
//				}
//
//			}
//		}
//		delete cell;
//	}
//
//	/*
//	 * store data results in a Logger object (see Logger.h)
//	 */
//
//	time_limit.stop();
//
//	Sat.copy(S);
//	Sat.set_data(covTime,covIt,covUnex,covSplit);
//	Sat.write_csv(pathname);
//	Sat.write_all(pathname, 100);
//
//	// print results
//	cout << "-- Stratégie : " << Sat.Strat << endl;
//	cout << "\t Temps de résolution : " << time_limit.get_time() << endl;
//	cout << "\t Nombre de solutions : " << S.size()<<endl<<endl;
//
//}

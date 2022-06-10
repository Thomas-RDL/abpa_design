#ifndef IBEX_TABLE2_CTC
#define IBEX_TABLE2_CTC

#include "ibex.h"

using namespace ibex;

/*
Table contractor for a minibex system
*/

class CtcTable: public Ctc{
public:
	CtcTable(const System& sys, const Matrix& M, const std::vector<std::string>& vars);
	CtcTable(const Matrix& M, const std::vector<int>& idx);

	CtcTable(const Matrix& M, const BitSet& bst);


	void contract(IntervalVector& box);

	virtual void contract(IntervalVector& box, ContractContext& context);


protected:
	std::vector<std::string> var_names; // variables names
	std::vector<int> var_indices; // map the intervalvector indices with the column indices of the matrix values
	const Matrix& values; // table of acceptable values
	std::vector<bool> cons; // consistante column values

};


inline void CtcTable::contract(IntervalVector& box){
	ContractContext context(box);
	contract(box,context);
}



inline CtcTable::CtcTable(const System& sys, const Matrix& M, const std::vector<std::string>& vars)
: Ctc(vars.size()), values(M){
	assert(M.nb_cols()<=sys.nb_var);
	assert(vars.size() == M.nb_cols());

	const std::vector<std::string> names = sys.var_names();
	// check how to implement with vector or matrix of variables
	// Search for var names in the system to map their index with M columns
	for (int i=0;i<vars.size();i++){
		for(int j=0;j<names.size();j++){
			if(names[j]==vars[i]){
				var_indices.push_back(j);
			}
		}
	}
}


inline CtcTable::CtcTable(const Matrix& M, const std::vector<int>& idx)
: Ctc(idx.size()), values(M), var_indices(idx){
	assert(var_indices.size() == M.nb_rows());
}

inline CtcTable::CtcTable(const Matrix& M, const BitSet& bst)
: Ctc(bst.size()), values(M){
	for( auto bit = bst.begin(); bit != bst.end(); bit++){
		var_indices.push_back(bit.el);
	}
	assert(var_indices.size() == M.nb_rows());

}



inline void CtcTable::contract(IntervalVector& box, ContractContext& context){
	IntervalVector hull(box.size());
	cons.resize(values.nb_cols());

	// Check which column are constant
	for(int j=0; j<values.nb_cols();j++){
		cons[j] = true;
		for(int i=0; i<var_indices.size();i++){
			if((values[i][j] & box[var_indices[i]]).is_empty())cons[j] = false;
		}
	}

	// for each variable compute the interval hull containing the consistant values in the table
	for(int i=0; i<var_indices.size();i++){
		for(int j=0; j<values.nb_cols();j++){
			if(cons[j]){
				if(!hull[var_indices[i]].is_unbounded()){
					hull[var_indices[i]] |= values[i][j];
				}else{
					hull[var_indices[i]] = values[i][j];
				}
			}
		}
	}
//	cout <<"hull :"<<hull<<endl;
//	cout<<" box : "<<box<<endl;
	box &= hull;
//	cout<<" box : "<<box<<endl;

}

#endif

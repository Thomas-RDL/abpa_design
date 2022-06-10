#ifndef IBEX_INT_CTC
#define IBEX_INT_CTC

#include "ibex.h"

using namespace ibex;

class CtcIntdomain: public Ctc{
public:

	CtcIntdomain(const System& sys, const std::vector<std::string>& vars);

	CtcIntdomain(const std::vector<int>& idx);

	CtcIntdomain(const BitSet& bst);

	void contract(IntervalVector& box);

	virtual void contract(IntervalVector& box, ContractContext& context);


protected:
	std::vector<std::string> var_names; // variables names
	std::vector<int> var_indices; // map the intervalvector indices with the column indices of the matrix values

};


inline void CtcIntdomain::contract(IntervalVector& box){
	ContractContext context(box);
	contract(box,context);
}



inline CtcIntdomain::CtcIntdomain(const System& sys, const std::vector<std::string>& vars)
: Ctc(vars.size()){

	assert(sys.nb_var()>=vars.size());
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


inline CtcIntdomain::CtcIntdomain( const std::vector<int>& idx)
: Ctc(idx.size()),var_indices(idx){
}

inline CtcIntdomain::CtcIntdomain(const BitSet& bst)
: Ctc(bst.size()){
	for( auto bit = bst.begin(); bit != bst.end(); bit++){
		var_indices.push_back(bit.el);
	}

}



inline void CtcIntdomain::contract(IntervalVector& box, ContractContext& context){
	// for each variable compute the interval hull containing the consistant values in the table
	for(int i=0; i<var_indices.size();i++){
			if((box[var_indices[i]] &= integer(box[var_indices[i]] )).is_empty()){
				box.set_empty();
				context.output_flags.add(FIXPOINT);
				return;
			}
			else{
				box[var_indices[i]] = integer(box[var_indices[i]] );
			}
		}
}

#endif

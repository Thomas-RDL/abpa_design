//============================================================================
//                                  I B E X                                   
// File        : Firsteger First bisector
// Author      : Thomas Richard de Latour
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : July 2020
//============================================================================

#ifndef __IBEX_BSC_RATIO_H__
#define __IBEX_BSC_RATIO_H__

#include "ibex_Bsc.h"
#include "ibex_BitSet.h"
#include "ibex_RoundRobin.h"
#include "ibex_LargestFirst.h"
#include "ibex_System.h"
#include "ibex_SmearFunction.h"
#include "ibex_NoBisectableVariableException.h"

namespace ibex {

/**
 * \ingroup bisector
 *
 * \brief ratio list bisector for MINLP system
 *
 * See #choose_var(const Cell& cell) for details.
 *
 */
class BscRatio : public Bsc {
public:
	/**
	 * \brief Create a bisector with ratio First heuristic.
	 *
	 */
	BscRatio(const BitSet& is_first, Bsc& bsc_first, Bsc& bsc_last, double prec,int var_ratio, double ratio=Bsc::default_ratio());

	/**
	 * \brief Create a bisector with ratio First heuristic.
	 */
	BscRatio(const BitSet& is_first, Bsc& bsc_first, Bsc& bsc_last, const Vector& prec,int var_ratio, double ratio=Bsc::default_ratio());

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by Bsc::bisect(...)
	 *
	 */
	virtual BisectionPoint choose_var(const Cell& cell);

	/**
	 * \brief Ratio to choose the split point.
	 *
	 * Ratio between the diameters of the left and right parts of a bisected
	 * interval.
	 */
	const double ratio;

	int  var_ratio;
protected :

	Bsc& bsc_first;

	Bsc& bsc_last;

	BitSet is_first;

	int nb_first ;
};
/*============================================ inline implementation ============================================ */

inline BscRatio::BscRatio(const BitSet& is_first, Bsc& bsc_first, Bsc& bsc_last, double prec, int var_ratio, double ratio) : Bsc(prec), bsc_first(bsc_first), bsc_last(bsc_last), is_first(is_first), var_ratio(var_ratio), ratio(ratio) {
	nb_first =0;
}

inline BscRatio::BscRatio(const BitSet& is_first, Bsc& bsc_first, Bsc& bsc_last, const Vector& prec, int var_ratio , double ratio) : Bsc(prec), bsc_first(bsc_first), bsc_last(bsc_last), is_first(is_first), var_ratio(var_ratio), ratio(ratio) {
	nb_first =0;
}

inline BisectionPoint BscRatio::choose_var(const Cell& cell) {

	const IntervalVector& box=cell.box;

	BitSet is_last;//création du bitset des dernières variables à bissecter
	is_last.resize(box.size());
	IntervalVector Lastbox(box.size(),0);//boîte contenant uniquement les dernières variables à bissecter
	IntervalVector Firstbox(box.size(),0); //boîte contenant uniquement les premières variables à bissecter

	for(int i=0;i<box.size();i++){
		if(!is_first[i]){
			is_last.add(i);
			Lastbox[i]=box[i];
		}else{ Firstbox[i]=box[i];
		}
	}
	Cell Firstcell(Firstbox,is_first.next(cell.bisected_var),cell.depth);
	Cell Lastcell(Lastbox,is_last.next(cell.bisected_var),cell.depth);


	if(cell.bisected_var >=is_first.max())
		Firstcell.bisected_var = is_first.min();


	if(cell.bisected_var >=is_first.max())
		Lastcell.bisected_var = is_last.min();

	try {
		if((!is_first.empty()) && nb_first < var_ratio){
			nb_first +=1;
			return bsc_first.choose_var(Firstcell);

		}
		else{throw NoBisectableVariableException(); }
	}
	catch (NoBisectableVariableException&) {

		nb_first = 0;
		return bsc_last.choose_var(Lastcell);
	}


}

} // end namespace ibex

#endif // __IBEX_BSC_RATIO__H__

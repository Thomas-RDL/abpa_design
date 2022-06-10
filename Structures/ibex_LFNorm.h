//============================================================================
//                                  I B E X                                   
// File        : Largest First bisector
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 19, 2012
// Last Update : Dec 25, 2017
//============================================================================

#ifndef __IBEX_L_F_NORM_H__
#define __IBEX_L_F_NORM_H__

#include "ibex_Bsc.h"
#include "ibex_NoBisectableVariableException.h"

namespace ibex {

/**
 * \ingroup bisector
 *
 * \brief normalized largest-first bisector.
 *
 */
class LFNorm : public Bsc {
public:

	/**
	 * \brief Create a bisector with largest-first heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double). By default, 0 which means an endless uniform bisection process.
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
	LFNorm(double prec=0, double ratio=Bsc::default_ratio());

	/**
	 * \brief Create a bisector with largest first heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double).
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
	LFNorm(const Vector& prec, double ratio=Bsc::default_ratio());

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by Bsc::bisect(...)
	 */
	virtual BisectionPoint choose_var(const Cell& cell);

	/**
	 * \brief Ratio to choose the split point.
	 *
	 * Ratio between the diameters of the left and right parts of a bisected
	 * interval.
	 */
	const double ratio;
	
	void set_box(IntervalVector& box);

 protected :
	virtual bool nobisectable (const IntervalVector& box, int i) const ;
	
	IntervalVector* init_box ;
};

inline LFNorm::LFNorm(double prec, double ratio1) : Bsc(prec), ratio(ratio1) {
	init_box = NULL;
}

LFNorm::LFNorm(const Vector& prec, double ratio1) : Bsc(prec), ratio(ratio1) {
	init_box = NULL;
}

  bool LFNorm::nobisectable(const IntervalVector & box, int i) const {
    return too_small (box, i);
  }


inline BisectionPoint LFNorm::choose_var(const Cell& cell) {

	if(init_box == NULL)
		throw NoBisectableVariableException();

	const IntervalVector& box=cell.box;

	int var =-1;
	double l=0.0;
	for (int i=0; i< box.size(); i++)	{
	  if (!(nobisectable (box,i)) && ((*init_box)[i].diam() >= this->prec(i))){
	  		double box_norm = box[i].diam()/(*init_box)[i].diam();
			if (var==-1) {
				var=i;
				l = uniform_prec()? box_norm : (box_norm/prec(i));
			}
			else {
				double l_tmp = uniform_prec()? box_norm : (box_norm/prec(i));
				if (l_tmp>l) {
					var = i;
					l = l_tmp;
				}
			}
		}
	}

	if (var !=-1){
		return BisectionPoint(var,ratio,true);
	}
	else {
		throw NoBisectableVariableException();
	}

}

inline void LFNorm::set_box(IntervalVector& box){
	init_box = &box;
}

} // end namespace ibex

#endif // __IBEX_L_F_NORM_H__

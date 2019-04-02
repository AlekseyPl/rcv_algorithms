/*
 * WtsCalc.h
 *
 *  Created on: Nov 9, 2018
 *      Author: aplotnikov
 */

#ifndef INCLUDE_WTSCALC_H_
#define INCLUDE_WTSCALC_H_

#include "Types.h"
#include <algorithm>
#include <vector>
#include <memory>

class WtsCalc {
public:

	WtsCalc( const SignalParams& sp, bool debug = false);
	~WtsCalc( );


	void CalculateWts( const DMatrix& RefSigsR, const DMatrix& RefSigsI,
					   const DVector& SigR, const DVector& SigI, int k);
	void CalculatePrevStAndCumWts(DVector& NewCumWts, DVector& OldCumWts, IVector& PrevStsCurCol);

private:

	DataType type;
	bool debug;
	int M;
	int L;
	int log2M;

	int ML;
	int ML1;
	int ML2;

	int Nt;
	int RefSigShift;
	int BufRefSigShift;

	int NSyms;

	// Объявление вектора весов переходов (Wts - branch // weights) размером ML элементов
	DVector wts;

	void CalculateComplex(const DMatrix& RefSigsR, const DMatrix& RefSigsI,
						  const DVector& SigR, const DVector& SigI, int k);
	void CalculateReal(const DMatrix& RefSigsR, const DVector& SigR, int k);
};

#endif /* INCLUDE_WTSCALC_H_ */

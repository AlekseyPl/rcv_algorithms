/*
 * AlphaBetaCommon.h
 *
 *  Created on: Nov 1, 2017
 *      Author: aplotnikov
 */
#ifndef ALPHABETACOMMON_H_
#define ALPHABETACOMMON_H_
 
#include "Types.h"

class AlphaBetaCommon {
public:
	AlphaBetaCommon( unsigned N_, unsigned M_, unsigned L_, AlgorithmType type_);
	~AlphaBetaCommon();

protected:
	AlgorithmType		type; 

	unsigned			L;
	unsigned			N;
	unsigned			M;
	unsigned			log2M;

	unsigned			ML;
	unsigned			ML1;
	unsigned			ML2;
	unsigned			NL;
	unsigned			NL1;

    int        			SymbsPrState;
};

#endif

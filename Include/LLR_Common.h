/*
 * LLR_Common.h
 *
 *  Created on: Nov 8, 2017
 *      Author: aplotnikov
 */

#ifndef INCLUDE_LLR_COMMON_H_
#define INCLUDE_LLR_COMMON_H_

#include "Types.h"

class CommonLLR {
public:
	CommonLLR( int N, int M, int L, int NSymb, AlgorithmType type );
	~CommonLLR( );

	void	LLR2SymbsPr( double *LLR );
	void	LLR2SymbsLogPr( double *LLR );


	DMatrix& GetInSymbLLR( )
    {
        return SymbsPr;
    }

	DMatrix& GetInSymbLogLLR( )
    {
        return SymbsLogPr;
    }

protected:
	int			    N;
	int			    NSymb;
	int			    NNSymb;
	int			    M;
	int			    L;

	int			    ML;
	int			    ML1;
	int			    ML2;
	int			    NL;
	int			    NL1;
	int			    log2M;

	double 			PrZer;
	double			LLRZer;
	double			PrOne;
	double			LLROne;
 

	AlgorithmType	type;

	DMatrix         SymbsPr;
	DMatrix		    SymbsLogPr;
	DMatrix		    BitsPr;
	DMatrix		    BitsLogPr;
	IMatrix		    Map;
};



#endif /* INCLUDE_LLR_COMMON_H_ */

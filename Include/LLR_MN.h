/*
 * LLR_MN.h
 *
 *  Created on: Nov 8, 2017
 *      Author: aplotnikov
 */

#ifndef INCLUDE_LLR_MN_H_
#define INCLUDE_LLR_MN_H_

#include "LLR_Common.h"

class MnLLR : public CommonLLR {
public:
	MnLLR( int N, int M, int L, int NSymb, AlgorithmType type, int MN_ );
	~MnLLR( );

    void    CalculateBitsLLR( double *OutLLR, StMatrix Alpha, StMatrix Beta );
	void    CalculateBitsMaxLLR(double *OutLLR, StMatrix Alpha, StMatrix Beta) ;

	StMatrix&	GetLambda()
	{
		return StLambda;
	}

private:

	int				MN;
    StMatrix        StLambda;
};



#endif /* INCLUDE_LLR_MN_H_ */

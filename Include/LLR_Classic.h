#ifndef LLR_CLASSIC_H
#define LLR_CLASSIC_H

#include "LLR_Common.h"

class ClassicLLR : public CommonLLR {
public:
	ClassicLLR( int N, int M, int L, int NSymb, AlgorithmType type );
	~ClassicLLR( );

	void	CalculateBitsLLR( double *OutLLR, DMatrix Alpha, DMatrix Beta );
	void	CalculateBitsMaxLLR( double *OutLLR, DMatrix Alpha, DMatrix Beta );

private:

    DMatrix		    Lambda;
};

#endif

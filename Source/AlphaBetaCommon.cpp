#include "../Include/AlphaBetaCommon.h"

AlphaBetaCommon::AlphaBetaCommon(unsigned N_, unsigned M_, unsigned L_, AlgorithmType type_ ):
	L( L_ ), N( N_ ), M( M_ ), type( type_ )
{
	unsigned k = 2;
	log2M = 1;

	while (k < M) {
		k <<= 1;
		log2M++;
	}

	ML2 = 1;
	for ( int k = 0; k < L-2; ++k )
		ML2 <<= log2M;

	ML1 = ML2 << log2M;
	ML  = ML1 << log2M;
	NL = N + L;
	NL1 = NL - 1;

}

AlphaBetaCommon::~AlphaBetaCommon( )
{

}

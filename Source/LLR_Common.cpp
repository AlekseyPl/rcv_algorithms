/*
 * LLR_Common.cpp
 *
 *  Created on: Nov 8, 2017
 *      Author: aplotnikov
 */


#include "../Include/LLR_Common.h"
#include <math.h>

CommonLLR::CommonLLR( int N_, int M_, int L_, int NSymb_, AlgorithmType type_ ):
	N( N_ ), M( M_ ), L( L_ ), NSymb( NSymb_ ), type( type_ ),
	PrZer( 0 ), PrOne( 0 ), LLRZer( 0 ), LLROne( 0 )
{
	int k = 2;
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
	NNSymb = N * NSymb;

	Map.resize( log2M );
	for ( int k = 0; k < log2M; ++k )
		Map[ k ].resize( M );
 
	switch( type ) {
		case TrueBCJR:
		case TrueMBCJR:

			BitsPr.resize( log2M + log2M );
			for ( int k = 0; k < log2M+log2M; ++k)
				BitsPr[ k ].resize( NNSymb );

			SymbsPr.resize( M );
			for ( int k = 0; k < M; ++k ) {
				SymbsPr[ k ].resize( NNSymb );
				SymbsPr[ k ].assign( NNSymb, 1 );
			}
			break;
		case MaxBCJR:
		case MaxMBCJR:

			// Указатель на массив BitsLogPr
			BitsLogPr.resize( log2M + log2M );
			for ( int k = 0; k < log2M+log2M; ++k)
				BitsLogPr[ k ].resize( NNSymb );

			SymbsLogPr.resize( M );
			for ( int k = 0; k < M; ++k ) {
				SymbsLogPr[ k ].resize( NNSymb );
				SymbsLogPr[ k ].assign( NNSymb, 0 );
			}
			break;
	}
}

CommonLLR::~CommonLLR( )
{
	for ( int k = 0; k < log2M; ++k )
		Map[ k ].clear( );
	Map.clear( );

	switch( type ) {
		case TrueBCJR:
		case TrueMBCJR:
		for ( int k = 0; k < log2M+log2M; ++k)
			BitsPr[ k ].clear();
		BitsPr.clear( );

		for ( int k = 0; k < M; ++k )
			SymbsPr[ k ].clear( );
		SymbsPr.clear( );
		break;
	case MaxBCJR:
	case MaxMBCJR:
		for ( int k = 0; k < log2M+log2M; ++k )
			BitsLogPr[ k ].clear( );
		BitsLogPr.clear( );

		for ( int k = 0; k < M; ++k )
			SymbsLogPr[k].clear( );
		SymbsLogPr.clear( );
		break;
	}
}


///////////////////////////////////////////////////////////////////////////
void CommonLLR::LLR2SymbsPr( double *LLR )
{
// Преобразуем llr в вероятности появления 0 и 1 и
// расставим их в массиве BitsPr размером N*NSymb x 2*log2M
	for ( int k = 0; k < NNSymb; ++k ) {
		for ( int n = 0; n < log2M; ++n ){
			// Вероятность появления 1
			BitsPr[ n + log2M ][ k ] =	1 / (1 + exp( *LLR++ ));
			// Вероятность появления 0
			BitsPr[ n ][ k ] 		 = 1 -	BitsPr[ n + log2M ][ k ];
		}
	}

// Составим карту выбора элементов BitsPr для вычисления
// вероятностей символов. Размер карты M x log2M
	int Mask = 0;
	for ( int k = 0; k < M; ++k ){
		Mask = 1 << log2M;
		for ( int n = 0; n < log2M; ++n ){
			Mask >>= 1;
			Map[ n ][ k ] = n + ( n && Mask ) ? log2M : 0;
		}
	}

// Вычислим вероятности символов и запишем их в массив размером N*NSymb x M
	for ( int k = 0; k < NNSymb; ++k ) {
		for ( int n = 0; n < M; ++n ) {
			for ( int m = 0; m < log2M; ++m)
				SymbsPr[ n ][ k ] *= BitsPr[ Map[ m ][ n ] ][ k ];
    		}
	}
}

///////////////////////////////////////////////////////////////////////////
void CommonLLR::LLR2SymbsLogPr( double *LLR )
{
// Преобразуем llr в логарифмы вероятностей появления 0 и 1 и
// расставим их в массиве BitsLogPr размером N*NSymb x 2*log2M
	int	m = 0;
	for ( int k = 0; k < NNSymb; ++k) {
		for ( int n = 0; n < log2M; ++n ){
			// Логарифм вероятности появления 1 = -ln(e^llr + 1)
			if (*( LLR + m ) > 0) 
				BitsLogPr[ n + log2M ][ k ] = -*( LLR + m );
			else 
				BitsLogPr[ n + log2M ][ k ] = 0;
			// Логарифм вероятности появления 0
			BitsLogPr[ n ][ k ] = *(LLR + m) + 	BitsLogPr[ n + log2M ][ k ];
			// Инкрементируем индекс в массиве llr
			m++;
		}
	}
// Составим карту выбора элементов BitsLogPr для вычисления
// вероятностей символов. Размер карты M x log2M
	int Mask;
	for ( int k = 0; k < M; ++k ){
		Mask = 1 << log2M;
		for ( int n = 0; n < log2M; ++n ){
			Mask >>= 1;
			Map[ n ][ k ] = n + ( n && Mask ) ? log2M : 0;

		}
	}

// Вычислим логарифмы вероятностей символов и запишем их в массив размером N*NSymb x M
	for ( int k = 0; k < NNSymb; ++k ) {
		for ( int n = 0; n < M; ++n ) {
			for ( int m = 0; m < log2M; ++m )
				SymbsLogPr[ n ][ k ] +=	BitsLogPr[ Map[ m ][ n ] ][ k ];
		}
	}
}

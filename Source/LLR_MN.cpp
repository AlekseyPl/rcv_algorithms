/*
 * LLR_MN.cpp
 *
 *  Created on: Nov 8, 2017
 *      Author: aplotnikov
 */
#include "../Include/LLR_MN.h"
#include <math.h>
#include <algorithm>
MnLLR::MnLLR( int N_, int M_, int L_, int NSymb_, AlgorithmType type_, int MN_):
	CommonLLR( N_, M_, L_, NSymb_, type_ ),  MN( MN_ )
{
	if (MN > ML1)
		MN = ML1;

	StLambda.resize( NL ); 
	for ( auto& it : StLambda )
		it.resize( MN );
	
}

MnLLR::~MnLLR( )
{
    for ( int  k = 0; k < NL; ++k )
        StLambda[ k ].clear();
    StLambda.clear( );
}

void	MnLLR::CalculateBitsLLR(double* outLLR, StMatrix Alpha, StMatrix Beta)
{
	for( auto& i: StLambda )
		i.assign( MN, ValueState( 0, 0 ) );

	for ( int t = 1; t < N + 1; ++t ) {
		for ( int k = 0; k < MN; ++k ) {
			for ( int i = 0; i < MN; ++i ) {
				if ( Alpha[ t ][ k ].state == Beta[ t ][ i ].state) {
					StLambda[ t ][ k ].value = Alpha[ t ][ k ].value * Beta[ t ][ i ].value;
					StLambda[ t ][ k ].state = Alpha[ t ][ k ].state;
					break;
				}
			}
		}
	}	
	for ( int t = 1; t < N+1; ++t ) {
		for ( int m = 0; m < log2M; ++m ) {
			PrZer = 0; // По-разному в True Max MaxStar
			PrOne = 0; // По-разному в True Max MaxStar

            for( int k = 0; k < MN; ++k ) {
                if( StLambda[ t ][ k ].state < ML2 )    PrZer += StLambda[ t ][ k ].value;
                else                                    PrOne += StLambda[ t ][ k ].value;
            }
			*outLLR++  = log(PrZer / PrOne); 
		}
	}
}

void	MnLLR::CalculateBitsMaxLLR( double* outLLR, StMatrix Alpha, StMatrix Beta )
{
	for( auto& i: StLambda )
		i.assign( MN, ValueState( 0, MInf) );

	for (int t = 1; t < N + 1; ++t) {
		for (int k = 0; k < MN; ++k) {
			for (int i = 0; i < MN; ++i) {
				if ( Alpha[ t ][ k ].state == Beta[ t ][ i ].state ) {
					StLambda[ t ][ k ].value = Alpha[ t ][ k ].value + Beta[ t ][ i ].value;
					StLambda[ t ][ k ].state = Alpha[ t ][ k ].state;
					break;
				}
			}
		}
	}

	for ( int t = 1; t < N+1; ++t ) {
		for ( int m = 0; m < log2M; ++m ) {
			LLRZer = MInf; // По-разному в True Max MaxStar
			LLROne = MInf; // По-разному в True Max MaxStar

            for( int k = 0; k < MN; ++k ) {
				if (StLambda[ t ][ k ].state < ML2)		LLRZer = std::max(LLRZer, StLambda[ t ][ k ].value);
				else 									LLROne = std::max(LLROne, StLambda[ t ][ k ].value);
            }
		 
			*outLLR++  = LLRZer - LLROne;
		}
	}
}

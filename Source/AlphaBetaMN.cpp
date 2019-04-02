/*
 * AlphaBetaMN.cpp
 *
 *  Created on: Nov 1, 2017
 *      Author: aplotnikov
 */
#include "../Include/AlphaBetaMN.h"
#include <algorithm>
#include <math.h>
#include <string.h>

AlphaBetaMN::AlphaBetaMN( int N_, int M_, int L_, AlgorithmType type_, int MN_, int Nt_, double Var_ ):
	AlphaBetaCommon( N_,M_,L_,type_ ), MN( MN_ ), Nt( Nt_ ), Var( Var_ )
{
	if (MN > ML1)  MN = ML1;

	Alpha.resize( NL );
	for ( int k = 0; k < NL; ++k )
		Alpha[k].resize( MN );
	 
	Beta.resize( NL );
	for ( int k = 0; k < NL; ++k )
		Beta[k].resize( MN );


	if( type == TrueMBCJR ) {
		Alpha[ 0 ][ 0 ].value  = 1;
		Beta[ NL1 ][ 0 ].value = 1;
	}
	else if ( type == MaxMBCJR ) {
		Alpha[ 0 ][ 0 ].value  = 0;
		Beta[ NL1  ][ 0 ].value = 0;
	}

	StateBuf.resize( M * MN );
 }

AlphaBetaMN::~AlphaBetaMN()
{
	for (int k = 0; k < NL; ++k)
		Alpha[ k ].clear();
	for (int k = 0; k < NL; ++k)
		Beta[ k ].clear();
	Alpha.clear( );
	Beta.clear( );
	
	StateBuf.clear( );
}

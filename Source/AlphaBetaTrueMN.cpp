/*
 * AlphaBetaTrueMN.cpp
 *
 *  Created on: Nov 15, 2017
 *      Author: aplotnikov
 */
#include "../Include/AlphaBetaTrueMN.h"
#include <math.h>
#include <string.h>
#include <algorithm>

AlphaBetaTrueMN::AlphaBetaTrueMN(int N_, int M_, int L_, AlgorithmType type_, int MN_, int Nt_, double Var_ ):
	AlphaBetaMN( N_, M_, L_, type_, MN_, Nt_, Var_ )
{

}
 
 
///////////////////////////////////////////////////////////////////////////
void AlphaBetaTrueMN::CalculateAlpha( DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift)
{
    double 	Sum;
	double 	R;
	double 	Buf;
	int 	PrevState;
	int		CurrState;
	SymbsPrState = symbPrShift;
	SpectrsPos	 = spectrsShift;
	RefShift = 0;
	MaxStateCount  	 	= 1;   // Расчёт альфы начинается с нулевого состояния
	ValueState tmp( 0, 0 );
 
	for ( int t = 1; t < L; ++t ){
		std::fill( Alpha[ t ].begin(), Alpha[ t ].end(), tmp);
		for ( int st = 0; st < MaxStateCount; ++st ) {
			for ( int m = 0; m < M; ++m ) { 
				PrevState = Alpha[ t - 1 ][ st ].state;
				CurrState = BackwardStates[ m ][ PrevState ];
				int RefOffset = PrevState;
				if (CurrState >= ML2) RefOffset += ML1;
				R = 0;
				for( int n = 0; n < Nt; ++n ) { 
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				}

				StateBuf[ st * M + m ].value = Alpha[ t - 1 ][ st ].value * exp(R * Var) * InSymbPr[ m ][ SymbsPrState ] ;
				StateBuf[ st * M + m ].state = CurrState;
			}
		}
		SetBestStartMPaths( Alpha[ t ], MaxStateCount );

		MaxStateCount <<= log2M;
		if (MaxStateCount > MN) MaxStateCount = MN;
		RefShift += ML;
		SpectrsPos++;
		SymbsPrState++;
	}

	for ( int t = L; t <= N; ++t ) {
		std::fill(Alpha[ t ].begin(), Alpha[ t ].end(), tmp);
		for ( int st = 0; st < MaxStateCount; ++st ) {
			for ( int m = 0; m < M; ++m ) { 
				PrevState = Alpha[ t - 1 ][ st ].state;
				CurrState = BackwardStates[ m ][ PrevState ];
				int RefOffset = PrevState;
				if (CurrState >= ML2) RefOffset += ML1;
				R = 0;
				for( int n = 0; n < Nt; ++n ) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				}

				StateBuf[ st * M + m ].value = Alpha[ t - 1 ][ st ].value * exp( R * Var ) * InSymbPr[ m ][ SymbsPrState ];
				StateBuf[ st * M + m ].state = CurrState;
			}
		}
		SetBestMPaths( Alpha[ t ], MaxStateCount );
		SpectrsPos++;
		SymbsPrState++;
	}

    int EndStateCount = ML2;
	for ( int t = N + 1; t < NL; ++t ) {
		std::fill(Alpha[ t ].begin(), Alpha[ t ].end(), tmp);
        MaxStateCount = ( EndStateCount > MN ) ? MN : EndStateCount;
		EndStateCount >>= log2M;
		RefShift += ML;

		for ( int st = 0; st < MaxStateCount; ++st ) {
			PrevState = Alpha[ t - 1 ][ st ].state;
	 		CurrState = BackwardStates[ 0 ][ PrevState ];
			int RefOffset = PrevState;
			if (CurrState >= ML2) RefOffset += ML1;
			R = 0;
			for (int n = 0; n < Nt; ++n) {
				Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
				R += Buf * Buf;
			}

			StateBuf[ st ].value = Alpha[ t - 1 ][ st ].value * exp( R * Var );
			StateBuf[ st ].state = CurrState;
		}
		SpectrsPos++;
		SetBestEndMPaths( Alpha[ t ], MaxStateCount );
	}

}

void AlphaBetaTrueMN::CalculateBeta( DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift)
{
	double R;
	double Buf;
	int RefOffset;
	MaxStateCount 	= 1;
	ValueState tmp(0, 0);
	int CurrState =	0; 
	int NextState = 0;

	RefShift = ML * ( 2 * L - 2 );
	SpectrsPos = spectrsShift + N + L - 2;

	for ( int t = N+L-2; t >= N; --t ) {
		for ( int st = 0; st < MaxStateCount; ++st ){
			for (int m = 0; m < M; ++m) { 
				NextState = Beta[ t + 1 ][ st ].state;
				CurrState =	ForwardStates[ m ][ NextState ];
				RefOffset = CurrState;
				if (NextState >= ML2)  RefOffset += ML1;
				R = 0;
				for (int n = 0; n < Nt; ++n) { 
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				} 

				StateBuf[ st * M + m ].value = Beta[ t + 1 ][ st ].value * exp( R * Var );
				StateBuf[ st * M + m ].state = CurrState;
			}
		}

		SetBestStartMPaths( Beta[ t ], MaxStateCount );
 		MaxStateCount <<= log2M;
 		if( MaxStateCount > MN ) MaxStateCount = MN;
 		RefShift 	  -= ML;
		SpectrsPos--;
	}
	SymbsPrState = symbPrShift + N - 1;
	for ( int t = N-1; t >= L-1; --t ) {
		fill(Beta[ t ].begin(), Beta[ t ].end(), tmp );
		for( int st = 0; st < MaxStateCount; ++st ) {
			for( int m = 0; m < M; ++m ) {
				NextState = Beta[ t + 1 ][ st ].state;
				CurrState = ForwardStates[ m ][ NextState ];
				RefOffset = CurrState;
				if (NextState >= ML2)  RefOffset += ML1;
				R = 0;
				for (int n = 0; n < Nt; ++n ) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				} 

				StateBuf[ st * M + m ].value = Beta[ t + 1 ][ st ].value * exp( R * Var ) * InSymbPr[ m ][ SymbsPrState ];
				StateBuf[ st * M + m ].state = CurrState;
			}
		}
		SpectrsPos--;
		SymbsPrState--;
		SetBestMPaths( Beta[ t ], MaxStateCount );
 	}

 	int EndStateCount = ML2;
	for ( int t = L - 2; t >= 0; --t ){
		
		fill(Beta[ t ].begin(), Beta[ t ].end(), tmp);
		RefShift  -= ML;
		MaxStateCount = ( EndStateCount > MN ) ? MN : EndStateCount;
		for ( int st = 0; st < MaxStateCount;  ++st ){
			NextState = Beta[ t + 1 ][ st ].state;
			CurrState = ForwardStates[ 0 ][ NextState ];
			RefOffset = CurrState;
			if (NextState >= ML2)  RefOffset += ML1;
			R = 0;
			for (int n = 0; n < Nt; ++n) {
				Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
				R += Buf * Buf;
			} 

			StateBuf[ st ].value = Beta[ t + 1 ][ st ].value * exp( R * Var ) * InSymbPr[ 0 ][ SymbsPrState ];
			StateBuf[ st ].state = CurrState;
		}
		SetBestEndMPaths( Beta[ t ], MaxStateCount);
		EndStateCount >>= log2M;
		SymbsPrState--;
		SpectrsPos--;
	}
}

void AlphaBetaTrueMN::SetBestStartMPaths(StVector& stVector, int count)
{
	std::sort(StateBuf.begin(), StateBuf.begin() + M * count, [ ](const ValueState& a, const ValueState& b) { return a.value < b.value; });

	int nextCount = ( count < MN ) ? ( count << 1 ) : MN;
	int k = M * count - 1;
	for (int i = 0; i < nextCount; ++i) 
		stVector[ i ] = StateBuf[ k-- ];

	Normalize( stVector, nextCount);
}

void AlphaBetaTrueMN::SetBestMPaths( StVector& stVector, int count )
{
    for( int i = 0; i < M * MN; ++i ) {
        for( int j = 0; j < M * MN; ++j ) {
			if (i != j) {
				if (StateBuf[ i ].state == StateBuf[ j ].state) {
					StateBuf[ i ].value += StateBuf[ j ].value;
					StateBuf[ j ].value = 0;
				}
			}
        }
    }
 
    std::sort( StateBuf.begin( ), StateBuf.begin() + M * MN, [] ( const ValueState& a, const ValueState& b ) { return a.value < b.value; } );

	int k = MN * M - 1;
	for (int i = 0; i < MaxStateCount; ++i)
		stVector[ i ] = StateBuf[ k-- ];

	Normalize( stVector, MaxStateCount );
}
void AlphaBetaTrueMN::SetBestEndMPaths(StVector& stVector, int count)
{
	for (int i = 0; i < count; ++i) {
		for (int j = 0; j < count; ++j) {
			if (i != j) {
				if (StateBuf[ i ].state == StateBuf[ j ].state) {
					StateBuf[ i ].value += StateBuf[ j ].value;
					StateBuf[ j ].value = 0;
				}
			}
		}
	}

	std::sort(StateBuf.begin(), StateBuf.begin() + count, [ ](const ValueState& a, const ValueState& b) { return a.value < b.value; });
	  
	int k = count - 1;
	for (int i = 0; i < count; ++i)
		stVector[ i ] = StateBuf[ k-- ];

	Normalize( stVector, count);
}


void AlphaBetaTrueMN::Normalize( StVector& stVector, int count )
{
	double	Sum        = 0;
	std::for_each(stVector.begin( ), stVector.begin( ) + count, [ &Sum ]( ValueState Beta ){ Sum += Beta.value; } );
	if ( Sum )
		std::for_each(stVector.begin( ), stVector.begin( ) + count, [ Sum ]( ValueState &Beta ){ Beta.value /= Sum; } );

}


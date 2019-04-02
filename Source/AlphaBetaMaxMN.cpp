/*
 * AlphaBetaMaxpMN.cpp
 *
 *  Created on: Nov 15, 2017
 *      Author: aplotnikov
 */

#include "../Include/AlphaBetaMaxMN.h"
#include <math.h>
#include <algorithm>
#include <string.h>
 

AlphaBetaMaxMN::AlphaBetaMaxMN(int N_, int M_, int L_, AlgorithmType type_, int MN_, int Nt_, double Var_) :
	AlphaBetaMN(N_, M_, L_, type_, MN_, Nt_, Var_)
{
 
}
///////////////////////////////////////////////////////////////////////////
void AlphaBetaMaxMN::CalculateAlpha(DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift)
{
	double 			Sum;
	double 			R;
	double 			Buf;
	unsigned		RefOffset;
	unsigned 		PrevState;
	unsigned		CurrState;
	SymbsPrState  = symbPrShift;
	SpectrsPos 	  = spectrsShift;
	RefShift 	  = 0;
	MaxStateCount = 1;   // Расчёт альфы начинается с нулевого состояния

	for (int t = 1; t < L; ++t) {
		for (int st = 0; st < MaxStateCount; ++st) {
			for (int m = 0; m < M; ++m) {
				PrevState = Alpha[ t - 1 ][ st ].state;
				CurrState = GetCurrStateFromPrev( PrevState, m );
				RefOffset = GetRefOffsetFromPrev( CurrState, PrevState );
 				R = 0;
				for (int n = 0; n < Nt; ++n) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				}
				
				StateBuf[ st * M + m ].value = Alpha[ t - 1 ][ st ].value + R * Var + InSymbPr[ m ][ SymbsPrState ];
				StateBuf[ st * M + m ].state = CurrState;
			}
		}
		SetBestStartMPaths(Alpha[ t ], MaxStateCount);

		MaxStateCount <<= log2M;
		if (MaxStateCount > MN) MaxStateCount = MN;
		RefShift += ML;
		SpectrsPos++;
		SymbsPrState++;
	}

	for (int t = L; t <= N; ++t) {
		for (int st = 0; st < MN; ++st) {
			for (int m = 0; m < M; ++m) {
				PrevState = Alpha[ t - 1 ][ st ].state;
				CurrState = GetCurrStateFromPrev( PrevState, m );
				RefOffset = GetRefOffsetFromPrev(CurrState, PrevState);
				R = 0;
				for (int n = 0; n < Nt; ++n) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				} 

				StateBuf[ st * M + m ].value = Alpha[ t - 1 ][ st ].value + R * Var + InSymbPr[ m ][ SymbsPrState ];
				StateBuf[ st * M + m ].state = CurrState;
			}
		}
		SetBestMPaths(Alpha[ t ], MaxStateCount);
		SpectrsPos++;
		SymbsPrState++;
	}

	int EndStateCount = ML2;
	for (int t = N + 1; t < NL; ++t) {

		fill(Alpha[ t ].begin(), Alpha[ t ].end(), ValueState(0, MInf));
		MaxStateCount = ( EndStateCount > MN ) ? MN : EndStateCount;
		EndStateCount >>= log2M;
		RefShift += ML;
		for (int st = 0; st < MaxStateCount; ++st) {
			PrevState = Alpha[ t - 1 ][ st ].state;
			CurrState = GetCurrStateFromPrev(PrevState, 0);
			RefOffset = GetRefOffsetFromPrev(CurrState, PrevState);
			R = 0;
			for (int n = 0; n < Nt; ++n) {
				Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
				R += Buf * Buf;
			} 

			StateBuf[ st ].value = Alpha[ t - 1 ][ st ].value + R * Var;
			StateBuf[ st ].state = CurrState;
		}
		SpectrsPos++;
		SetBestEndMPaths(Alpha[ t ], MaxStateCount);
	}
}

void AlphaBetaMaxMN::CalculateBeta(DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift)
{
	double		R;
	double		Buf;
	unsigned	RefOffset;
	unsigned	CurrState;
	unsigned	NextState;
 
	MaxStateCount = 1;
	RefShift	  = ML * ( 2 * L - 2 );
	SpectrsPos	  = spectrsShift + N + L - 2;

	for (int t = N + L - 2; t >= N; --t) {
		for (int st = 0; st < MaxStateCount; ++st) {
			for (int m = 0; m < M; ++m) {
				NextState = Beta[ t + 1 ][ st ].state;
				CurrState = GetCurrStateFromNext( NextState, m );
				RefOffset = GetRefOffsetFromNext( CurrState, NextState );
 
				R = 0;
				for (int n = 0; n < Nt; ++n) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				} 

				StateBuf[ st * M + m ].value = Beta[ t + 1 ][ st ].value + R * Var;
				StateBuf[ st * M + m ].state = CurrState;
			}
		}

		SetBestStartMPaths(Beta[ t ], MaxStateCount);
		MaxStateCount <<= log2M;
		if (MaxStateCount > MN) MaxStateCount = MN;
		RefShift -= ML;
		SpectrsPos--;
	}
	SymbsPrState = symbPrShift + N - 1;
	for (int t = N - 1; t >= L - 1; --t) { 
		for (int st = 0; st < MaxStateCount; ++st) {
			for (int m = 0; m < M; ++m) {
				NextState = Beta[ t + 1 ][ st ].state;
				CurrState = GetCurrStateFromNext( NextState, m );
				RefOffset = GetRefOffsetFromNext(CurrState, NextState);

				R = 0;
				for (int n = 0; n < Nt; ++n) {
					Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
					R += Buf * Buf;
				} 

				StateBuf[ M * st + m ].value = Beta[ t + 1 ][ st ].value + R * Var + InSymbPr[ m ][ SymbsPrState ];
				StateBuf[ M * st + m ].state = CurrState;
			}
		}
		SpectrsPos--;
		SymbsPrState--;
		SetBestMPaths(Beta[ t ], MaxStateCount);
	}

	int EndStateCount = ML2;
	for (int t = L - 2; t >= 0; --t) {
		fill(Beta[ t ].begin(), Beta[ t ].end(), ValueState(0,MInf));
		RefShift -= ML;
		MaxStateCount = ( EndStateCount > MN ) ? MN : EndStateCount;
		for (int st = 0; st < MaxStateCount; ++st) {

			NextState = Beta[ t + 1 ][ st ].state;
			CurrState = GetCurrStateFromNext( NextState, 0 );
			RefOffset = GetRefOffsetFromNext(CurrState, NextState);

			R = 0;
			for (int n = 0; n < Nt; ++n) {
				Buf = *( RefSigs[ RefShift + RefOffset ] + n ) - *( Spectrs[ SpectrsPos ] + n );
				R += Buf * Buf;
			}

			StateBuf[ st ].value = Beta[ t + 1 ][ st ].value + R * Var + InSymbPr[ 0 ][ SymbsPrState ];
			StateBuf[ st ].state = CurrState;
		}
		SetBestEndMPaths(Beta[ t ], MaxStateCount);
		EndStateCount >>= log2M;
		SymbsPrState--;
		SpectrsPos--;
	}
}
 
void AlphaBetaMaxMN::SetBestStartMPaths(StVector& stVector, int count)
{
	std::sort(StateBuf.begin(), StateBuf.begin() + M * count, [ ](const ValueState& a, const ValueState& b) { return a.value < b.value; });

	int nextCount = ( count < MN ) ? ( count << 1 ) : MN;
	int k = M * count - 1;
	for (int i = 0; i < nextCount; ++i)
		stVector[ i ] = StateBuf[ k-- ];

	Normalize( stVector, nextCount );
}

void AlphaBetaMaxMN::SetBestMPaths(StVector& stVector, int count)
{
	for (int i = 0; i < M * MN; ++i) {
		for (int j = 0; j < i; ++j) {
			if (StateBuf[ i ].state == StateBuf[ j ].state) {
				StateBuf[ i ].value = std::max(StateBuf[ i ].value, StateBuf[ j ].value);
				StateBuf[ j ].value = MInf;
			}
		}
		for( int j = i + 1; j < M * MaxStateCount; ++j ) {
			if (StateBuf[ i ].state == StateBuf[ j ].state) {
				StateBuf[ i ].value = std::max(StateBuf[ i ].value, StateBuf[ j ].value);
				StateBuf[ j ].value = MInf;
			}
		}
	}

	std::sort(StateBuf.begin(), StateBuf.begin() + M * MN, [ ](const ValueState& a, const ValueState& b) { return a.value > b.value; });
	memcpy(&stVector[ 0 ], &StateBuf[ 0 ], MN * sizeof( ValueState ));
 
	Normalize(stVector, MN);
}

void AlphaBetaMaxMN::SetBestEndMPaths(StVector& stVector, int count)
{
	for (int i = 0; i < count; ++i) {
		for (int j = 0; j < count; ++j) {
			if (i != j) {
				if (StateBuf[ i ].state == StateBuf[ j ].state) {
					StateBuf[ i ].value = std::max(StateBuf[ i ].value, StateBuf[ j ].value);
					StateBuf[ j ].value = MInf;
				}
			}
		}
	}

	std::sort(StateBuf.begin(), StateBuf.begin() + count, [ ](const ValueState& a, const ValueState& b) { return a.value > b.value; });
	memcpy(&stVector[ 0 ], &StateBuf[ 0 ], count * sizeof(ValueState));

 

	Normalize( stVector, count);
}

void AlphaBetaMaxMN::Normalize(StVector& stVector, int count)
{
	double	Max = MInf;
	std::for_each( stVector.begin(), stVector.begin() + count, [ &Max ]( ValueState vs ) { if( Max < vs.value ) Max = vs.value; } );
	std::for_each( stVector.begin(), stVector.begin() + count, [ Max ]( ValueState &vs ) { vs.value -= Max; } );
}




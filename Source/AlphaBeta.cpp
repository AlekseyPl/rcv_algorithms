#include "../Include/AlphaBeta.h"
#include <math.h>
#include <algorithm>
#include <string.h>

AlphaBeta::AlphaBeta( int N_, int M_,  int L_, int MN_, AlgorithmType type_ ):
	AlphaBetaCommon( N_,M_,L_,type_ )
{
	double Buf = 0;
	MInf = -1/Buf;


	// Р�РЅРёС†РёР°Р»РёР·Р°С†РёСЏ СѓРєР°Р·Р°С‚РµР»РµР№ РЅР° РјР°СЃСЃРёРІС‹ Alpha, Beta Рё Lambda
    // (СЂР°Р·РјРµСЂ ML1 x NL), РІРєР»СЋС‡Р°СЏ РёРЅРёС†РёР°Р»РёР·Р°С†РёСЋ Р·РЅР°С‡РµРЅРёР№ Alpha,
    Alpha.resize( NL );
	for ( int k = 0; k < NL; ++k )
        Alpha[k].resize( ML1 );

    Beta.resize( NL );
	for ( int k = 0; k < NL; ++k )
		Alpha[k].resize( ML1 );

 
}

AlphaBeta::~AlphaBeta( )
{
 
    for ( int  k = 0; k < NL; ++k )
        Alpha[ k ].clear( );
	for ( int k = 0; k < NL; ++k )
        Beta[ k ].clear( );
    Alpha.clear( );
    Beta.clear( );

}


void AlphaBeta::ReInit( )
{

	if( type == TrueBCJR ) {
	    for( int i = 0; i < NL; ++i ) {
	    	Alpha[ i ].assign( ML1, 0 );
	    	Beta[ i ].assign(  ML1, 0 );
	    }
		Alpha[ 0 ][ 0 ]  = 1;
		Beta[ NL1 ][ 0 ] = 1;
    }
	else if ( type == MaxBCJR ) {
	    for( int i = 0; i < NL; ++i ) {
	    	Alpha[ i ].assign( ML1, MInf );
	    	Beta[ i ].assign(  ML1, MInf );
	    }
		Alpha[ 0 ][ 0 ]  = 0;
		Beta[ NL1 ][ 0 ] = 0;
	}
}

void AlphaBeta::CalculateAlpha( DMatrix R, DMatrix SymbsPr, int SymbsPrShift )
{
	if( type == TrueBCJR )
		CalcTrueAlpha( R, SymbsPr, SymbsPrShift );
	else if( type == MaxBCJR )
		CalcMaxAlpha( R, SymbsPr, SymbsPrShift );

}

void AlphaBeta::CalculateBeta( DMatrix R, DMatrix SymbsPr, int SymbsPrShift )
{
	if( type == TrueBCJR )
		CalcTrueBeta( R, SymbsPr, SymbsPrShift );
	else if( type == MaxBCJR )
		CalcMaxBeta( R, SymbsPr, SymbsPrShift );
}


void AlphaBeta::CalcTrueAlpha( DMatrix R, DMatrix SymbsPr, int SymbsPrShift )
{
 
	int CurState, PrevState, RefSigNum, SymbNum;
	int CurStateStep, PrevStateStep, RefSigNumStep;
	int SymbNumCounter, SymbNumCounterMax;
	int PrevState0, MaxCurState;
 
	double Sum;

	CurStateStep	  = ML2;
	PrevStateStep     = ML1;
	RefSigNumStep	  = PrevStateStep;
	SymbNumCounterMax = 1;
	SymbsPrState	  = SymbsPrShift;
	for ( int t = 1; t < L; ++t ){
		PrevState      = 0;
		RefSigNum      = 0;
		SymbNum        = 0;
		SymbNumCounter = 0;
		Sum            = 0; 
		for ( int CurState = 0; CurState < ML1; CurState += CurStateStep ){

			Alpha[ t ][ CurState ] = Alpha[ t - 1 ][ PrevState ] * R[ t - 1 ][ RefSigNum ] * SymbsPr[ SymbNum ][ SymbsPrState ];
	
			PrevState = ( PrevState + PrevStateStep ) % ML1;
			RefSigNum += RefSigNumStep;

			if ( ++SymbNumCounter == SymbNumCounterMax) {
				SymbNum++;
				SymbNumCounter = 0;
			}
		}
		std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + ML1, [ &Sum ]( double alpha ){ Sum += alpha; } );
		std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + ML1, [ Sum ]( double& alpha ){ alpha /= Sum; } );

		PrevStateStep = CurStateStep;
		RefSigNumStep = PrevStateStep;
		CurStateStep >>= log2M;
		SymbNumCounterMax <<= log2M;
		SymbsPrState++;
	}

	for ( int t = L; t <= N; ++t ){
		PrevState0     = 0;
		RefSigNum      = 0;
		SymbNum        = 0;
		SymbNumCounter = 0;
		Sum            = 0; // РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
		for ( DVector::iterator it = Alpha[ t ].begin(); it < Alpha[ t ].begin() + ML1; ++it ){
			for ( int n = 0; n < M; ++n )
				*it += Alpha[ t - 1 ][ PrevState0 + n ] * R[ t - 1 ][ RefSigNum++ ] * SymbsPr[ SymbNum ][ SymbsPrState ];


			PrevState0 = ( PrevState0 + M ) % ML1;

			if ( ++SymbNumCounter == ML2 ){
				SymbNum++;
				SymbNumCounter = 0;
			}
		}

		std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + ML1, [ &Sum ]( double alpha ){ Sum += alpha; } );
		if( Sum ) std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + ML1, [ Sum ]( double& alpha ){ alpha /= Sum; } );

		SymbsPrState++;
	}

	MaxCurState = ML2;
	for ( int t = N + 1; t < NL; ++t ){
		PrevState0 = 0;
		RefSigNum  = 0;
		Sum        = 0; // пїЅпїЅ-пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ True Max MaxStar
		for ( int CurState = 0; CurState < MaxCurState; ++CurState ){
			for( int n = 0; n < M; ++n ) 
				Alpha[ t ][ CurState ] += Alpha[ t - 1 ][ PrevState0 + n ] * R[ t - 1 ][ RefSigNum++ ];
		
			PrevState0 += M;
		}
		std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + MaxCurState, [ &Sum ]( double alpha ){ Sum += alpha; } );
		std::for_each( Alpha[ t ].begin(), Alpha[ t ].begin() + MaxCurState, [ Sum ]( double& alpha ){ alpha /= Sum; } );
		MaxCurState >>= log2M;
	}
}
 
///////////////////////////////////////////////////////////////////////////
void AlphaBeta::CalcTrueBeta( DMatrix R, DMatrix SymbsPr, int SymbsPrShift )
{
	int SymbsPrState;
	int CurState, NextState, NextState0, RefSigNum, RefSigNum0;
	int MaxCurState, NextState0Step, RefSigNum0Step;
	int NextStateCount, NextState0Count;
	 
	double Sum; 
	MaxCurState = M;
	for ( int t = N+L-2; t >= N; --t ) {
		NextState      = 0;
		NextStateCount = 0;
		RefSigNum      = 0;
		Sum            = 0; 
		for ( int CurState = 0; CurState < MaxCurState; ++CurState ) {
			// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStarr
			Beta[ t ][ CurState ] = Beta[ t + 1 ][ NextState ] * R[ t ][ RefSigNum++ ];
			// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
			if ( ++NextStateCount == M ) {
				NextStateCount = 0;
				NextState++;
			}
		}
		std::for_each( Beta[ t ].begin(), Beta[ t ].begin() + MaxCurState, [ &Sum ]( double Beta ){ Sum += Beta; } );
		std::for_each( Beta[ t ].begin(), Beta[ t ].begin() + MaxCurState, [ Sum ]( double &Beta ){ Beta /= Sum; } );
		MaxCurState <<= log2M;
	}

	NextState0Count = 0;
	SymbsPrState = SymbsPrShift + N-1;
	for ( int t = N-1; t >= L-1; --t ) {
		NextState0 = 0;
		RefSigNum0 = 0;
		Sum        = 0; // РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
		for ( int CurState = 0; CurState < ML1; ++CurState ) {
			NextState = NextState0;
			RefSigNum = RefSigNum0++;
			for ( int n = 0; n < M; ++n ) {
				// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
				Beta[ t ][ CurState ] += Beta[ t + 1 ][ NextState ] * R[ t ][ RefSigNum ] * SymbsPr[ n ][ SymbsPrState ];

				// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
				NextState += ML2;
				RefSigNum += ML1;
			}
			// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
			if ( ++NextState0Count == M ) {
				NextState0Count = 0;
				NextState0++;
			}
		}
		std::for_each( Beta[ t ].begin(), Beta[ t ].begin() + ML1, [ &Sum ]( double Beta ){ Sum += Beta; } );
		if( Sum ) std::for_each( Beta[ t ].begin(), Beta[ t ].begin() + ML1, [ Sum ]( double &Beta ){ Beta /= Sum; } );
		SymbsPrState--;
	}

	NextState0Step = 1;
	RefSigNum0Step = M;
	for ( int t = L - 2; t >= 0; --t ){
		NextState0 = 0;
		RefSigNum0 = 0;
		Sum        = 0; // РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
		for ( int CurState = 0; CurState < ML1; CurState += RefSigNum0Step ){
			NextState = NextState0;
			RefSigNum = RefSigNum0;
		 	for ( int n = 0; n < M; ++n ) {
				// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
				Beta[ t ][ CurState ] += Beta[ t + 1 ][ NextState ] * R[ t ][ RefSigNum ] * SymbsPr[ n ][ SymbsPrState ];
			
				// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
				NextState += ML2;
				RefSigNum += ML1;
			}
			// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
			Sum += Beta[ t ][ CurState ];

			// РџРѕ-СЂР°Р·РЅРѕРјСѓ РІ True Max MaxStar
			RefSigNum0 += RefSigNum0Step;
			NextState0 += NextState0Step;
		}

		for ( int CurState = 0; CurState < ML1; CurState += RefSigNum0Step )
			Beta[ t ][ CurState ] /= Sum; // пїЅпїЅ-пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ True Max MaxStar

		NextState0Step = RefSigNum0Step;
		RefSigNum0Step <<= log2M;
		SymbsPrState--;
	}
}

void AlphaBeta::CalcMaxAlpha( DMatrix R, DMatrix SymbsPr, int SymbsPrShift )
{
    int SymbsPrState;
    int CurState, PrevState, RefSigNum, SymbNum;
    int CurStateStep, PrevStateStep, RefSigNumStep;
    int SymbNumCounter, SymbNumCounterMax;
    int PrevState0, MaxCurState;

    // Начало разной реализации в True Max MaxStar
    double Buf, Max;

    Buf = 0;

    CurStateStep  = ML2;
    PrevStateStep = ML1;
    RefSigNumStep = PrevStateStep;
    SymbNumCounterMax = 1;
    SymbsPrState = SymbsPrShift;
    for ( int t = 1; t < L; ++t ){
        PrevState      = 0;
        RefSigNum      = 0;
        SymbNum        = 0;
        SymbNumCounter = 0;
        Max            = MInf; // По-разному в True Max MaxStar
        for ( int CurState = 0; CurState < ML1; CurState += CurStateStep){

        	Alpha[ t ][ CurState ] = Alpha[ t - 1][ PrevState ] + R[ t - 1 ][ RefSigNum ] + SymbsPr[ SymbNum ][ SymbsPrState ];
            if ( Alpha[ t ][ CurState ] > Max)
                Max = Alpha[ t ][ CurState ];

            PrevState = ( PrevState + PrevStateStep ) % ML1;
            RefSigNum += RefSigNumStep;

            if ( ++SymbNumCounter == SymbNumCounterMax){
                SymbNum++;
                SymbNumCounter = 0;
            }
        }

        for ( int CurState = 0; CurState < ML1; CurState += CurStateStep)
            Alpha[ t ][ CurState ] -= Max; // По-разному в True Max MaxStar

        PrevStateStep = CurStateStep;
        RefSigNumStep = PrevStateStep;
        CurStateStep >>= log2M;
        SymbNumCounterMax <<= log2M;
        SymbsPrState++;
    }

    for ( int t = L; t <= N; ++t ){
        PrevState0     = 0;
        RefSigNum      = 0;
        SymbNum        = 0;
        SymbNumCounter = 0;
        Max            = MInf; // По-разному в True Max MaxStar
        for ( int CurState = 0; CurState < ML1; ++CurState ){
        	for ( int n = 0; n < M; ++n ){
                // Начало разной реализации в True Max MaxStar
                Buf = Alpha[ t - 1][ PrevState0 + n ] + R[ t - 1 ][ RefSigNum++ ] + SymbsPr[ SymbNum ][ SymbsPrState ];
                if (Buf > Alpha[ t ][ CurState ] )
					Alpha[ t ][ CurState ] = Buf;
                // Конец разной реализации в True Max MaxStar
            }

			PrevState0 = ( PrevState0 + M ) % ML1;
            if ( ++SymbNumCounter == ML2){
                SymbNum++;
                SymbNumCounter = 0;
            }
        }

        std::for_each( Alpha[ t ].begin( ), Alpha[ t ].end( ), [ &Sum ]( double alpha ){ Max = ( alpha > Max) ? alpha : Max ;})
        std::for_each( Alpha[ t ].begin( ), Alpha[ t ].end( ), [ Sum ]( double& alpha ){ alpha -= Max; } )

        SymbsPrState++;
    }

    MaxCurState = ML2;
    for ( int t = N+1; t < N+L; ++t ){
        PrevState0 = 0;
        RefSigNum  = 0;
        Max        = MInf; // По-разному в True Max MaxStar
        for ( int CurState = 0; CurState < MaxCurState; ++CurState ){
            for ( int n = 0; n < M; ++n ){
                // Начало разной реализации в True Max MaxStar
                Buf = Alpha[ t - 1 ][ PrevState0 + n ] + R[ t - 1 ][ RefSigNum++ ];
                if (Buf > Alpha[ t ][ CurState ] )
                	Alpha[ t ][ CurState ] = Buf;
                // Конец разной реализации в True Max MaxStar
            }
             PrevState0 += M;
        }

        std::for_each( Alpha[ t ].begin( ), Alpha[ t ].begin( ) + MaxCurState, [ &Sum ]( double alpha ){ Max = ( alpha > Max) ? alpha : Max ;})
		std::for_each( Alpha[ t ].begin( ), Alpha[ t ].begin( ) + MaxCurState, [ Sum ]( double& alpha ){ alpha -= Max; } )


        MaxCurState >>= log2M;
    }
}

///////////////////////////////////////////////////////////////////////////
void AlphaBeta::CalcMaxBeta( double **R,  double **SymbsPr, int SymbsPrShift )
{
    int SymbsPrState;
    int CurState, NextState, NextState0, RefSigNum, RefSigNum0;
    int MaxCurState, NextState0Step, RefSigNum0Step;
    int NextStateCount, NextState0Count;

    // Начало разной реализации в True Max MaxStar
    double Buf, Max;
    double MInf;

    Buf = 0;
    MInf = -1/Buf;
    // Конец разной реализации в True Max MaxStar

    MaxCurState = M;
    for ( int t = N+L-2; t >= N; --t ) {
        NextState      = 0;
        NextStateCount = 0;
        RefSigNum      = 0;
        Max            = MInf; // По-разному в True Max MaxStar
        for ( int CurState = 0; CurState < MaxCurState; ++CurState ) {
            // Начало разной реализации в True Max MaxStar
            Beta[ t ][ CurState ] = Beta[ t + 1 ][ NextState ] + R[ t ][ RefSigNum++ ];
             // Конец разной реализации в True Max MaxStar

            if ( ++NextStateCount == M ) {
                NextStateCount = 0;
                NextState++;
            }
        }

        std::for_each( Beta[ t ].begin( ), Beta[ t ].end( ), [ &Sum ]( double beta ){ Max = ( beta > Max) ? alpha : Max ;})
		std::for_each( Beta[ t ].begin( ), Beta[ t ].begin( ) + MaxCurState, [ Sum ]( double& beta ){ beta -= Max; } )

        MaxCurState <<= log2M;
    }

    NextState0Count = 0;
    SymbsPrState = SymbsPrShift + N-1;
    for ( int t = N-1; t >= L-1; --t ) {
        NextState0 = 0;
        RefSigNum0 = 0;
        Max        = MInf; // По-разному в True Max MaxStar
        for ( int CurState = 0; CurState < ML1; ++CurState ) {
            NextState = NextState0;
            RefSigNum = RefSigNum0;
            for ( int n = 0; n < M; ++n ) {
                // Начало разной реализации в True Max MaxStar
                Buf = Beta[ t + 1 ][ NextState ] + R[ t ][ RefSigNum ] + SymbsPr[ n ][ SymbsPrState ];
                if (Buf > Beta[ t ][ CurState ] )
                	Beta[ t ][ CurState ] = Buf;
                // Конец разной реализации в True Max MaxStar
                NextState += ML2;
                RefSigNum += ML1;
            }
            if ( ++NextState0Count == M) {
                NextState0Count = 0;
                NextState0++;
            }
            RefSigNum0++;
        }
        std::for_each( Beta[ t ].begin( ), Beta[ t ].end( ), [ &Sum ]( double beta ){ Max = ( beta > Max) ? alpha : Max ;})
		std::for_each( Beta[ t ].begin( ), Beta[ t ].end( ), [ Sum ]( double& beta ){ beta -= Max; } )

        SymbsPrState--;
    }

    NextState0Step = 1;
    RefSigNum0Step = M;
    for ( int t = L-2; t >= 0; --t){
        NextState0 = 0;
        RefSigNum0 = 0;
        Max        = MInf; // По-разному в True Max MaxStar
        for ( DVector::iterator it = 0; it < ML1; it += RefSigNum0Step ){
            NextState = NextState0;
            RefSigNum = RefSigNum0;

            for ( int n = 0; n < M; ++n ) {
                // Начало разной реализации в True Max MaxStar
                Buf = Beta[ t + 1 ][ NextState ] + R[ t ][ RefSigNum ] + SymbsPr[ n ][ SymbsPrState ];
                if (Buf > *it )
                	*it = Buf;
                // Конец разной реализации в True Max MaxStar
                NextState += ML2;
                RefSigNum += ML1;
            }
            // Начало разной реализации в True Max MaxStar
            if ( *it > Max)
                Max = *it;
            // Конец разной реализации в True Max MaxStar
            RefSigNum0 += RefSigNum0Step;
            NextState0 += NextState0Step;
        }

        for (CurState = 0; CurState < ML1; CurState += RefSigNum0Step)
        	Beta[ t ][ CurState ] -= Max; // По-разному в True Max MaxStar

        NextState0Step = RefSigNum0Step;
        RefSigNum0Step <<= log2M;
        SymbsPrState--;
    }

}



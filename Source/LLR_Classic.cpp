#include "../Include/LLR_Classic.h"
#include <math.h>


ClassicLLR::ClassicLLR( int N_, int M_, int L_, int NSymb_, AlgorithmType type_ ):
	CommonLLR( N_, M_, L_, NSymb_, type_ )
{
    Lambda.resize( NL );
	for ( int k = 0; k < NL; ++k )
		Lambda[ k ].resize( ML1 );
}

ClassicLLR::~ClassicLLR( )
{
    for ( int  k = 0; k < NL; ++k )
    	Lambda[ k ].clear();
    Lambda.clear( );
}
///////////////////////////////////////////////////////////////////////////
void ClassicLLR::CalculateBitsLLR( double *OutLLR, DMatrix Alpha, DMatrix Beta  )
{
 	int NumStates, StateStep, StateShift;
	 
	double PrZer, PrOne;

	for ( int t = 1; t < N + 1; ++t )
		for ( int k = 0; k < ML1; ++k )
			Lambda[ t ][ k ] = Alpha[ t ][ k ] * Beta[ t ][ k ]; // ��-������� � True Max MaxStar

	StateStep = ML1 >> log2M;
	for ( int t = 1; t < L-1; ++t ) {
		NumStates = (ML1 >> 1);
		for ( int m = 0; m < log2M; ++m ) {
			PrZer = 0; // По-разному в True Max MaxStar
			PrOne = 0; // По-разному в True Max MaxStar
			StateShift = 0;
			while ( StateShift < ML1 ){
				for ( int k = StateShift; k < StateShift + NumStates; k += StateStep) {
					PrZer += Lambda[ t ][ k ];				// По разному в True Max MaxStar
					PrOne += Lambda[ t ][ k + NumStates ];	// По-разному в True Max MaxStar
				}
				StateShift += 2 * NumStates;
			}
			*OutLLR++ = ( ( PrOne > 0.0 ) ) ? log(PrZer / PrOne) : 0.0; // ��-������� � True Max MaxStar
			NumStates >>= 1;
		}
		StateStep >>= log2M;
	}
	
	for ( int t = L-1; t < N+1; ++t ) {
		NumStates = (ML1 >> 1);
		for ( int m = 0; m < log2M; ++m ) {
			PrZer = 0; // По-разному в True Max MaxStar
			PrOne = 0; // По-разному в True Max MaxStar
			StateShift = 0;
			while ( StateShift < ML1 ) {
				for ( int k = StateShift; k < StateShift + NumStates; ++k ) {
					PrZer += Lambda[ t ][ k ];				// По разному в True Max MaxStar
					PrOne += Lambda[ t ][ k + NumStates ];	// По-разному в True Max MaxStar
				}
				StateShift += 2 * NumStates;
			}
			*OutLLR++  = ( ( PrOne > 0.0 ) ) ? log(PrZer / PrOne) : 0.0; // ��-������� � True Max MaxStar
			NumStates >>= 1;
		}
	}

}

///////////////////////////////////////////////////////////////////////////
void ClassicLLR::CalculateBitsMaxLLR( double *OutLLR,  DMatrix Alpha, DMatrix Beta )
{
	int NumStates, StateStep, StateShift;
	// ������ ������ ���������� � True Max MaxStar
	double LLRZer, LLROne;

	LLRZer = 0;
	double MInf = -1/LLRZer;
	// ����� ������ ���������� � True Max MaxStar

	for ( int t = 1; t < N + 1; ++t )
		for ( int k = 0; k < ML1; ++k )
			Lambda[ t ][ k ] = Alpha[ t ][ k ] + Beta[ t ][ k ]; // ��-������� � True Max MaxStar

	StateStep = ML1 >> log2M;
	for ( int t = 1; t < L-1; ++t ) {
		NumStates = (ML1 >> 1);
		for (int m = 0; m < log2M; ++m ) {
			LLRZer = MInf; // ��-������� � True Max MaxStar
			LLROne = MInf; // ��-������� � True Max MaxStar
			StateShift = 0;
			while (StateShift < ML1){
				for ( int k = StateShift; k < StateShift + NumStates; k += StateStep) {
					if ( Lambda[ t ][ k ] > LLRZer) // ��-������� � True Max MaxStar
						LLRZer = Lambda[ t ][ k ]; // ��-������� � True Max MaxStar

					if ( Lambda[ t ][ k + NumStates ] > LLROne) // ��-������� � True Max MaxStar
						LLROne = Lambda[ t ][ k + NumStates ]; // ��-������� � True Max MaxStar
				}
				StateShift += 2 * NumStates;
			}
			*OutLLR++ = LLRZer - LLROne; // ��-������� � True Max MaxStar
			NumStates >>= 1;
		}
		StateStep >>= log2M;
	}

	for ( int t = L - 1; t < N + 1; ++t ) {
		NumStates = (ML1 >> 1);
		for ( int m = 0; m < log2M; ++m) {
			LLRZer = MInf; // ��-������� � True Max MaxStar
			LLROne = MInf; // ��-������� � True Max MaxStar
			StateShift = 0;
			while (StateShift < ML1){
				for ( int k = StateShift; k < StateShift + NumStates; k++) {
					if ( Lambda[ t ][ k ] > LLRZer) // ��-������� � True Max MaxStar
						LLRZer = Lambda[ t ][ k ]; // ��-������� � True Max MaxStar

					if ( Lambda[ t ][ k + NumStates ] > LLROne) // ��-������� � True Max MaxStar
						LLROne = Lambda[ t ][ k + NumStates ]; // ��-������� � True Max MaxStar
				}
				StateShift += 2 * NumStates;
			}
			*OutLLR++ = LLRZer - LLROne; // ��-������� � True Max MaxStar
			NumStates >>= 1;
		}
	}
}




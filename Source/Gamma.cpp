#include "../Include/Gamma.h"
#include <math.h>
 
Gamma::Gamma( int N_, int M_, int L_, int Nt_, double Var_, bool isSigComplex_ ):
	N( N_ ), M( M_ ), L( L_ ), Nt( Nt_ ), Var( Var_ ), isSigComplex( isSigComplex_ ),
	ML1( 1 )
{
	int k = 2;
	log2M = 1;

	while ( k < M ) {
		k <<= 1;
		log2M++;
	}

	for ( int k = 0; k < L-1; ++k )
		ML1 <<= log2M; // ML1 = ML2 * M;
	// ML = M^L
	ML  = ML1 << log2M; // ML = ML1 * M;
	// NL = N + L
	NL = N + L;
	// NL1 = N + L - 1
	NL1 = NL - 1;

    R.resize( NL1 );
	for ( auto it = R.begin(); it < R.end(); ++it )
		it->resize( ML );
}

Gamma::~Gamma( )
{
	for ( auto it = R.begin(); it < R.end(); ++it )
		it->clear( );
	R.clear();
}

void Gamma::ReInit( )
{
	for( auto it = R.begin(); it < R.end( ); ++it )
		it->assign( it->size( ), 0 );
}
///////////////////////////////////////////////////////////////////////////
void Gamma::CalculateTrueR( double **SpectrsR, double **SpectrsI, double **RefSigsR,	double **RefSigsI, int SpectrsShift )
// Функция вычисления экспонент от отношений эвклидовых расстояний к
// удвоенной дисперсии, взятых со знаком минус
{ 
	PreCalculateR( SpectrsR, RefSigsR, SpectrsShift );
	if ( isSigComplex )
		PreCalculateR( SpectrsI, RefSigsI, SpectrsShift );

 	int MStep = ML1;
	for ( int k = 0; k < L-1;  ++k ) {
		for ( int m = 0; m < ML; m += MStep)
			R[ k ][ m ] = exp( R[ k ][ m ] * Var);
		MStep >>= log2M;
	}

	for ( int k = L-1; k < N; ++k ) {
		for ( int m = 0; m < ML; ++m )
			R[ k ][ m ] = exp( R[ k ][ m ] * Var);
	}

	int NL1 = N + L - 1;
	int MMax = ML1;
	for ( int k = N; k < NL1; ++k ) {
		for ( int m = 0; m < MMax; ++m )
			R[ k ][ m ] = exp( R[ k ][ m ] * Var);
		MMax >>= log2M;
	}
}

void Gamma::CalculateMaxR( double **SpectrsR, double **SpectrsI, double **RefSigsR, double **RefSigsI, int SpectrsShift )
// Функция вычисления логарифма экспонент от отношений эвклидовых расстояний к
// удвоенной дисперсии, взятых со знаком минус
{
    int MStep, MMax;

	PreCalculateR( SpectrsR, RefSigsR, SpectrsShift );
	if ( isSigComplex )
		PreCalculateR( SpectrsI, RefSigsI, SpectrsShift );
     

    MStep = ML1;
	for ( int k = 0; k < L-1;  ++k ) {
		for ( int m = 0; m < ML; m += MStep)
			R[ k ][ m ] = ( R[ k ][ m ] * Var);
        MStep >>= log2M;
    }

	for ( int k = L-1; k < N; ++k ) {
		for ( int m = 0; m < ML; ++m )
			R[ k ][ m ] = ( R[ k ][ m ] * Var);
    }

    MMax = ML1;
    for ( int k = N; k < N+L-1; ++k) {
        for ( int m = 0; m < MMax; ++m)
        	R[ k ][ m ] = ( R[ k ][ m ] * Var);
        MMax >>= log2M;
    }
 
}

///////////////////////////////////////////////////////////////////////////
// Расчёт евклидовых расстояний символов принятого сигнала и опорных символов
void Gamma::PreCalculateR( double **Spectrs, double **RefSigs, int SpectrsShift  )
{
	int RefShift = 0;
	int SpectrsPos = SpectrsShift;
	int MStep = ML1;
	double Buf = 0;

    // Расчёт для первых L-1 ТИ
	for ( int k = 0; k < L-1; ++k ) {
		for ( int m = 0; m < ML; m += MStep) {
			for ( int n = 0; n < Nt; ++n ) {
				Buf = *( RefSigs[ RefShift + m ] + n ) -
						*( Spectrs[ SpectrsPos ] + n );
				R[ k ][ m ] += Buf * Buf;
			}
		}
		SpectrsPos++;
		MStep >>= log2M;
		RefShift += ML;
	}

	 // Расчёт для средних N-L+1 ТИ
	for ( int k = L-1; k < N; ++k ) {
		for ( int m = 0; m < ML; ++m ) {
			for ( int n = 0; n < Nt; ++n ) {
				Buf = *( RefSigs[ RefShift + m ] + n ) -
						*( Spectrs[ SpectrsPos ] + n );
				R[ k ][ m ] += Buf * Buf;
			}
		}
		SpectrsPos++;
	}

	 // Расчёт для последних L-1 ТИ
	RefShift += ML;
	int MMax = ML1;
	int NL1 = N + L - 1;
	for ( int k = N; k < NL1;  ++k ) {
		for ( int m = 0; m < MMax; ++m ) {
			for ( int n = 0; n < Nt; ++n ) {
				Buf = *( RefSigs[ RefShift + m ] + n ) -
						*( Spectrs[ SpectrsPos ] + n );
				R[ k ][ m ] += Buf * Buf;
			}
		}
		SpectrsPos++;
		MMax >>= log2M;
		RefShift += ML;
	}
}

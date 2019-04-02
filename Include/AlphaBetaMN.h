/*
 * AlphaBetaMN.h
 *
 *  Created on: Nov 1, 2017
 *      Author: aplotnikov
 */
#ifndef ALPHABETAMN_H_
#define ALPHABETAMN_H_


#include "AlphaBetaCommon.h"

class AlphaBetaMN : public AlphaBetaCommon {
public:

	AlphaBetaMN( int N, int M, int L, AlgorithmType type, int MN, int Nt, double Var );
	virtual ~AlphaBetaMN( );

	virtual	void	CalculateAlpha( DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift) = 0;
	virtual	void	CalculateBeta(  DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift) = 0;

    StMatrix&    GetAlpha( ) 
    {
        return Alpha;
    }
    StMatrix&    GetBeta( ) 
    {
        return Beta;
    }
	 
protected:
	unsigned			MN;

//	for gamma calc
	unsigned			Nt;
	double				Var;

	StMatrix			Alpha;
	StMatrix			Beta;

	StVector			StateBuf;

	unsigned			RefShift;
	unsigned			MaxStateCount;
	unsigned			SpectrsPos;
 

    virtual	void		Normalize( StVector& stVector, int count ) = 0;
	virtual	void		SetBestStartMPaths(StVector& stVector, int count) = 0;
	virtual	void		SetBestMPaths(StVector& stVector, int count) = 0;
	virtual	void		SetBestEndMPaths(StVector& stVector, int count) = 0;

	inline	unsigned	GetCurrStateFromPrev( unsigned PrevState, unsigned dec )
	{
		return  ( PrevState >> 1 ) + dec * ML2;
	}
	inline	unsigned	GetCurrStateFromNext( unsigned NextState, unsigned dec )
	{
		return	M * NextState % ML1 + dec;
	}

	inline unsigned		GetRefOffsetFromNext(unsigned CurrState, unsigned NextState)
	{
		unsigned	RefOffset = CurrState;
		if (NextState >= ML2)  RefOffset += ML1;
		return RefOffset;
	}

	inline unsigned		GetRefOffsetFromPrev(unsigned CurrState, unsigned PrevState)
	{
		unsigned	RefOffset = PrevState;
		if (CurrState >= ML2) RefOffset += ML1;
		return RefOffset;
	}
};


#endif 

/*
 * AlphaBetaTrueMN.h
 *
 *  Created on: Nov 15, 2017
 *      Author: aplotnikov
 */

#ifndef INCLUDE_ALPHABETATRUEMN_H_
#define INCLUDE_ALPHABETATRUEMN_H_

#include "AlphaBetaMN.h"

class AlphaBetaTrueMN : public AlphaBetaMN  {
public:
	AlphaBetaTrueMN( int N, int M, int L, AlgorithmType type, int MN, int Nt, double Var );
	virtual ~AlphaBetaTrueMN() { };
	virtual	void	CalculateAlpha( DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift );
	virtual	void	CalculateBeta(  DMatrix InSymbPr, double **Spectrs, double **RefSigs, int symbPrShift, int spectrsShift );

private: 

	virtual	void		Normalize( StVector& stVector, int count );
	virtual	void		SetBestStartMPaths(StVector& stVector, int count);
	virtual	void		SetBestMPaths(StVector& stVector, int count);
	virtual	void		SetBestEndMPaths( StVector& stVector, int count);
};


#endif /* INCLUDE_ALPHABETATRUEMN_H_ */

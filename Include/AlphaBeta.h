#include "mex.h"
#include "AlphaBetaCommon.h"

class AlphaBeta : public AlphaBetaCommon {
public:

	AlphaBeta( int N, int M, int L, int MN, AlgorithmType type );
	~AlphaBeta( );

	void	ReInit( );
	
	void	CalculateAlpha( DMatrix R, DMatrix InSymbPr, int SymbsPrShift);
	void	CalculateBeta( DMatrix R,DMatrix InSymbPr, int SymbsPrShift );
 

	DMatrix& GetAlpha( )
	{
		return Alpha;
	}

	DMatrix& GetBeta( )
	{
		return Beta;
	}
private:

	DMatrix Alpha;
	DMatrix Beta;

	double	MInf;

	void	CalcMaxAlpha( DMatrix R, DMatrix InSymbPr, int SymbsPrShift );
	void	CalcTrueAlpha( DMatrix R, DMatrix InSymbPr, int SymbsPrShift );

	void	CalcMaxBeta( DMatrix R, DMatrix InSymbPr, int SymbsPrShift );
	void	CalcTrueBeta( DMatrix R, DMatrix InSymbPr, int SymbsPrShift );
};

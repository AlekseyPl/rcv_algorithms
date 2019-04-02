#include "Types.h"

class Gamma {
public:
	Gamma( int N, int M, int L, int Nt, double Var, bool isSigComplex );
	~Gamma( );

	void	ReInit( );

	void	CalculateTrueR( double **SpectrsR, double **SpectrsI, double **RefSigsR,	double **RefSigsI, int SpectrsShift );
    void	CalculateMaxR( double **SpectrsR, double **SpectrsI, double **RefSigsR,	double **RefSigsI, 	int SpectrsShift );

    DMatrix& GetR( )
    {
    	return R;
    }
private:

	
	int		Nt;
	int		N;
	int		M;
	int		L;

	int		ML;
	int		ML1;
	int		NL;
	int		NL1;
	int		log2M;
	bool	isSigComplex;
	double	Var;

	DMatrix R;

	void	PreCalculateR(  double **Spectrs, double **RefSigs, int SpectrsShift );
};

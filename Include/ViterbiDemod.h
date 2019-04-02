#ifndef VITERBIDEMOD_H
#define VITERBIDEMOD_H

#include "Types.h"
#include <memory>
 
class Sym2Bit;
class WtsCalc;

class ViterbiDemodulator {
public:
	ViterbiDemodulator( const SignalParams& sp, int TBLen, bool debug = false );
	~ViterbiDemodulator();

	void process(const DMatrix& RefSigsR,const DMatrix& RefSigsI, const DMatrix& SigR, const DMatrix& SigI,
				 DVector& OutLLR, int SigsShift, int BitNumShift);
private:

	bool debug;

	int TBLen;
	DataType type;
	ConstellationType constType;

	int N;
	int M;
	int L;
	int log2M;

	int ML;
	int ML1;
	int ML2;

	int Nt;
	int RefSigShift;
	int BufRefSigShift;

	int NSyms;
	int NBits;


	
// ��������� �� ������� ����� �����, �.�. ����������� ����� ���������,
// CumWts - cumulative (path) weights. NewCumWts � OldCumWts ����� ���
// ����, ����� �� ��������� ��������� ����� �����, �������� ��
// ���������� ���� ��������� ���� �������� ����� �� ������� ����.
// ��������� BufCumWts ����� ��� ������ �������� ���������� NewCumWts �
// OldCumWts ��� �������� � ������ ����.
	DVector NewCumWts;
	DVector OldCumWts;

// ������ ��� �������� �������� ������� ���������, �� ������� ������
// �������� ���� (previous states) (ML1 � TBLen).
	IMatrix PrevSts;

	std::shared_ptr< Sym2Bit > conv;
	std::shared_ptr< WtsCalc > wtsCalc;
};

#endif

#include "mex.h"
#include <memory>
#include <algorithm>
#include <vector>
#include <string.h>

#define VECTOR
#define VITERBI

#include "../Include/Types.h"
#include "../Include/ViterbiDemod.h"


using namespace std;

enum MexParamsIdx {
	L_idx		= 0,
	M_idx		= 1,
	N_idx		= 2,
	Nt_idx		= 3,
	Nsymb_idx   = 4,
	Var_idx     = 5, // not used in Viterbi
	Spectrs_idx = 6,
	RefSigs_idx = 7,
	InLLR_idx   = 8, // not used in Viterbi
	TBLen_idx   = 9
};

 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// �������-������� //////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // ��������! � ������� ��� ������� ������ �� ������, ������� ���
    // "���������" ������� ������ ����� ���� ��������� ���� � �������
    // ������� ����������!
    //
    // ������ ������ mex-������� �� MATLAB:
    //
    //   OutLLR = MexViterbi(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
    //          InLLR, TBLen);
    //      L     - ���-�� ���������, ����������� �����;
    //      M     - ������ ���������, ����������� �����;
    //      N     - ���-�� ������������� �������� � ����� NOFDM-�������,
    //              ����������� �����;
    //      Nt    - ���-�� �������� ������� �� ���� ��, ����������� �����;
    //      NSymb - ���-�� NOFDM-�������� � ����� �����, ����������� �����;
    //      Var   - ��������� ����, ������������ �����; - �� ������������
    //      Spectrs - (Nt*(N+L-1) x NSymb) ������ ������������ /
    //                ����������� �������� - ����� �� ����� ������������
    //                ��� ������ (Nt x (N+L-1)*NSymb);
    //      RefSigs - (Nt � (M^L)*(2*L-1)) ������ ������������ /
    //                ����������� ������� ��������;
    //      InLLR   - (N*NSymb*log2(M) � 1) ������ ������������ ��������
    //                ��������� llr ���; - �� ������������, ����� ��������
    //                ��, ��� ������;
    //      TBLen   - ��������� �������� �� 1 �� (N+L-1).
    //
    //   �������� ����������:
    //      OutLLR  - (N*NSymb*log2(M) � 1) ������ ������������ ��������
    //                ������������� llr ���.
    //
    // �������� � ��������� ����������! ��� �������� ���������� �������,
    // �.�. �������, �� Matlab � C++ ������ ������������ � ����������
    // ������, ������ ������� ���������� �� ������� - �� ��������!
    
    // ���������� ����������
	// ������� ����������
	// ������������� ����������-�������
		int NSymb;

	// ��������� �� ������������ (� ������) �����
		DMatrix SpectrsR;
		DMatrix SpectrsI;
		DMatrix RefSigsR;
		DMatrix RefSigsI;


	// �������� ����� ������ ������������, �.�. ���������� ��������
	// ���������� ������� �� ��������� ������� ���������� ��������
	// ������� ����������������� �������
		int TBLen;

// �������� ����������
		int k, SpectrsShift, BitNumShift;

		SignalParams params;
    // ��������� �������� ������� ����������
	// ���������� ���������
		params.L  = static_cast<int>(mxGetScalar(prhs[L_idx]));
	// ������ ���������
		params.M  = static_cast<int>(mxGetScalar(prhs[M_idx]));
	// ���-�� ������������� �������� � ����� NOFDM-�������
		params.N  = static_cast<int>(mxGetScalar(prhs[N_idx]));
	// ���-�� �������� ������� �� ���� ��
		params.Nt = static_cast<int>(mxGetScalar(prhs[Nt_idx]));
	// ���-�� NOFDM-�������� � ����� �����
		NSymb = static_cast<int>(mxGetScalar(prhs[Nsymb_idx]));
	// ��������� ���� - �� ������������!
	//    Var = (double)(mxGetScalar(prhs[5]));
	// TBLen
		TBLen = static_cast<int>(mxGetScalar(prhs[TBLen_idx]));

    // ������ ����� ������������ �������� - ��� ��������� ����    
         
		int log2M = 1;
		k = 2;
		while (k < params.M) {
			k <<= 1;
			log2M++;
		}
		int ML2 = 1;
		for( int k = 0; k < params.L-2; ++k )
			ML2 <<= log2M;

		params.log2M  = log2M;
		params.ML2    = ML2;
		params.ML1    = ML2 << params.log2M;
		params.ML     = params.ML1 << params.log2M;

		params.ML2L1  = params.ML * ( 2 * params.L - 1 );
		params.NL     = params.N + params.L;
		params.NL1    = params.NL -  1;
		params.NBits  = params.N * params.log2M;
	// NNSymb = N * NSymb
		int32_t NNSymb = params.N * NSymb;

	// ����, ����������� �� �� ����������� ������ ��� ���
		bool isSigComplex = mxIsComplex(prhs[Spectrs_idx]);
		params.type = isSigComplex ? DataType::complex : DataType::real;
		params.constType = ( params.M > 2 ) ? ConstellationType::m_ary : ConstellationType::binary;


// ����������� ��������� �������� ������� ����������

		SpectrsR.resize(params.NL1*NSymb);
		auto ptrR = mxGetPr(prhs[Spectrs_idx]);
		for( int i =0; i < SpectrsR.size(); ++i ) {
			SpectrsR.resize(params.Nt);
			for( auto& v : SpectrsR[ i ] ) v = *ptrR++;
		}


		if (isSigComplex) {	// ������ �����
			SpectrsI.resize(params.NL1*NSymb);
			auto ptrI = mxGetPi(prhs[Spectrs_idx]);
			for( int i =0; i < SpectrsI.size(); ++i ) {
				SpectrsI.resize(params.Nt);
				for( auto& v : SpectrsI[ i ] ) v = *ptrI++;
			}
		}

		RefSigsR.resize(params.ML2L1);
		auto ptrRefR = mxGetPr(prhs[RefSigs_idx]);
		for( int i =0; i < params.ML2L1; ++i ) {
			RefSigsR.resize(params.Nt);
			for( auto& v : RefSigsR[ i ] ) v = *ptrRefR++;
		}

		if (isSigComplex) {		// ������ �����
			RefSigsI.resize(params.ML2L1);
			auto ptrRefI = mxGetPi(prhs[RefSigs_idx]);
			for( int i =0; i < params.ML2L1; ++i ) {
				RefSigsI.resize(params.Nt);
				for( auto& v : RefSigsI[ i ] ) v = *ptrRefI++;
			}
		}


    // ��������� ������ ��� �������� ���������� - ������ OutLLR. ���
	// ��������� ������ ���� ������ ������ ��������, ����� MATLAB ����
    // �������� ���������.
		plhs[0] = mxCreateDoubleMatrix(NNSymb*params.log2M, 1, mxREAL);
		double *OutLLR = mxGetPr(plhs[0]);
		std::vector<double> vecOut;
		vecOut.reserve(NNSymb*params.log2M);

		shared_ptr< ViterbiDemodulator > vitD = make_shared< ViterbiDemodulator >( params, TBLen );

    // ���� �� ���������� NOFDM-��������
		SpectrsShift = 0; // ����� ������ �������� Spectrs �
		BitNumShift  = 0; // OutLLR
		for (k = 0; k < NSymb; k++) {
			// ����� �������� C-�������
				vitD->process( RefSigsR, RefSigsI, SpectrsR, SpectrsI, vecOut, SpectrsShift, BitNumShift);
			// ������� �������� ������ �� NOFDM-��������
				SpectrsShift += params.NL1; // ����� ������ �������� Spectrs �
				BitNumShift  += params.N*params.log2M; // OutLLR
		}
//		memcpy(OutLLR, &vecOut[0], vecOut.size() * sizeof(double));

//		for( auto&i : vecOut) 			*OutLLR++ = i;
//		for( int i = 0;)
        
    return;
}

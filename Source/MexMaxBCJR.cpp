#include "mex.h"
#include "../Include/Types.h"
#include "../Include/Gamma.h"
#include "../Include/AlphaBeta.h"
#include "../Include/LLR_Classic.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // ��������! � ������� ��� ������� ������ �� ������, ������� ���
    // "���������" ������� ������ ����� ���� ��������� ���� � �������
    // ������� ����������!
    //
    // ������ ������ mex-������� �� MATLAB:
    //
    //   OutLLR = MexMaxBCJR(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
    //          InLLR, TBLen);
    //
    //   ������� ����������:
    //      L     - ���-�� ���������, ����������� �����;
    //      M     - ������ ���������, ����������� �����;
    //      N     - ���-�� ������������� �������� � ����� NOFDM-�������,
    //              ����������� �����;
    //      Nt    - ���-�� �������� ������� �� ���� ��, ����������� �����;
    //      NSymb - ���-�� NOFDM-�������� � ����� �����, ����������� �����;
    //      Var   - ��������� ����, ������������ �����;
    //      Spectrs - (Nt*(N+L-1) x NSymb) ������ ������������ /
    //                ����������� �������� - ����� �� ����� ������������
    //                ��� ������ (Nt x (N+L-1)*NSymb);
    //      RefSigs - (Nt � (M^L)*(2*L-1)) ������ ������������ /
    //                ����������� ������� ��������;
    //      InLLR   - (N*NSymb*log2(M) � 1) ������ ������������ ��������
    //                ��������� llr ���;
    //      TBLen   - ��������� �������� �� 1 �� (N+L-1) - �� ������������!
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
                int L, M, N, Nt, NSymb;
                int TBLen;
            // ������������ ����������-�������
                double Var;
            // ��������� �� ������������ (� ������) �����(�)
            // ����������-��������
                double **SpectrsR, **SpectrsI, **RefSigsR, **RefSigsI,
                        *InLLR;

        // �������� ����������
            // ��������� �� ������ OutLLR
                double *OutLLR;

        // ���������� ����������
            // ����� ������������ ��������: log2M = log2(M), ML = M^L;
            // ML1 = M^(L-1); ML2 = M^(L-2); NL = N+L; NL1 = N+L-1,
            // NNSymb = N*NSymb
                int log2M, ML, ML1, ML2, NL, NL1, NNSymb;
            // ��������
                int k, SpectrsShift, InSymbsLogPrShift, BitNumberShift;
            // ����, ����������� �� �� ����������� ������ ��� ���
                bool isSigComplex;
            // R, Alpha, Beta, Lambda, InSymbsLogPr
                double **R, **Alpha, **Beta, **Lambda, **InSymbsLogPr;

    // ��������� �������� ������� ����������
        // ���������� ���������
            L  = (int)(mxGetScalar(prhs[0]));
        // ������ ���������
            M  = (int)(mxGetScalar(prhs[1]));
        // ���-�� ������������� �������� � ����� NOFDM-�������
            N  = (int)(mxGetScalar(prhs[2]));
        // ���-�� �������� ������� �� ���� ��
            Nt = (int)(mxGetScalar(prhs[3]));
        // ���-�� NOFDM-�������� � ����� �����
            NSymb = (int)(mxGetScalar(prhs[4]));
        // ��������� ����
            Var = (double)(mxGetScalar(prhs[5]));
        // TBLen - �� ������������!
        //    TBLen = (int)(mxGetScalar(prhs[9]));

    // ������ ����� ������������ �������� - ��� ��������� ����
        // log2M = log2(M)
            log2M = 1;
            k = 2;
            while (k < M) {
                k <<= 1;
                log2M++;
            }
        // ML2 = M^(L-2)
            ML2 = 1;
            for (k = 0; k < L-2; k++)
                ML2 <<= log2M; // ML2 *= M;
        // ML1 = M^(L-1)
            ML1 = ML2 << log2M; // ML1 = ML2 * M;
        // ML = M^L
            ML  = ML1 << log2M; // ML = ML1 * M;
        // NL = N + L
            NL = N + L;
        // NL1 = N + L - 1
            NL1 = NL - 1;
        // NNSymb = N * NSymb
            NNSymb = N * NSymb;
        // -1/(2*Var) - ��������� ���������� ����������
            Var = -1/(2*Var);

        // ����, ����������� �� �� ����������� ������ ��� ���
            isSigComplex = mxIsComplex(prhs[6]);
            
    // ����������� ��������� �������� ������� ����������
        // ��������� �� ������ Spectrs
            // ������������ �����
                SpectrsR = new double*[NL1*NSymb];
                SpectrsR[0] = mxGetPr(prhs[6]);
                for (k = 1; k < NL1*NSymb; k++)
                    SpectrsR[k] = SpectrsR[k-1] + Nt;
            // ������ �����
                if (isSigComplex) {
                    SpectrsI = new double*[NL1*NSymb];
                    SpectrsI[0] = mxGetPi(prhs[6]);
                    for (k = 1; k < NL1*NSymb; k++)
                        SpectrsI[k] = SpectrsI[k-1] + Nt;
                }
            
        // ��������� �� ������ RefSigs
            // ������������ �����
                RefSigsR = new double*[ML*(2*L-1)];
                RefSigsR[0] = mxGetPr(prhs[7]);
                for (k = 1; k < ML*(2*L-1); k++)
                    RefSigsR[k] = RefSigsR[k-1] + Nt;
            // ������ �����
                if (isSigComplex) {                
                    RefSigsI = new double*[ML*(2*L-1)];
                    RefSigsI[0] = mxGetPi(prhs[7]);
                    for (k = 1; k < ML*(2*L-1); k++)
                        RefSigsI[k] = RefSigsI[k-1] + Nt;
                }

// ��������� �� ������ InSymbsPr
	InLLR = mxGetPr(prhs[8]);


	AlgorithmType   aType = TrueMBCJR;
	MType           mType = bestBeta;
	//MType mType = ( MType )MNtype;

	AlphaBeta* alphabeta = new	AlphaBeta( N, M, L, aType );
	ClassicLLR*	 cllr 	 = new	ClassicLLR( N, M, L, NSymb, aType );
	Gamma*		 gamma   = new 	Gamma( N, M, L, Nt, Var, isSigComplex  )

    // ��������� ������ ��� �������� ���������� - ������ OutLLR. ���
	// ��������� ������ ���� ������ ������ ��������, ����� MATLAB ����
    // �������� ���������.
        plhs[0] = mxCreateDoubleMatrix(NNSymb*log2M, 1, mxREAL);
        OutLLR = mxGetPr(plhs[0]);
            
    // ��������� ������ ��� ���������� R, Alpha, Beta, Lambda, InSymbsLogPr
    // ���� �� ����� ��������� ���� ��������, �� ��� ��� ����������
    // ����������, ����� ��� ����������� ��� �������� ����������
    // ��� ��������� ������ ����� double ���������� ���������� ������ �����
    // ���� �����, ��� ��������� ����� mxCreateDoubleMatrix ������������
    // ���������� ������! �������� �������� ���, ��� R ���� �������� ������
    // ��� ������� �������; Alpha, Beta - ����� ���������������� ������
    // ���� �������� � ������ ���� ���; Lambda - �� ����� ����������������;
    // InSymbsLogPr - �������� ���������������� ������ LLR2SymbsLogPr


    // ������ ���� ��������� - �������������� ������� llr � ���������
    // ������������ �������� ��� ����� �����
        cllr->LLR2SymbsLogPr( InLLR );
        DMatrix&	SymbLogPr = cllr->GetInSymbLogLLR( );

        SpectrsShift      = 0; // ����� ������ �������� Spectrs,
        InSymbsLogPrShift = 0; // InSymbsLogPr �
        BitNumberShift    = 0; // OutLLR
        for (k = 0; k < NSymb; ++k ) {     // ���� �� ���������� NOFDM-��������
            // ������ R - �������� ������������ �������� ������� �������� �
            // �������� �������
        		gamma->ReInit( );
        		gamma->CalculateMaxR( Spectrs, SpectrsI, RefSigR, RefSigI, SpectrsShift );
        		DMatrix& R = gamma->GetR( );

            // ������ ������� Alpha
                alphabeta->CalculateAlpha( R, SymbLogPr, SymbsPrShift );

            // ������ ������� Beta
                alphabeta->CalculateBeta( R, SymbLogPr, SymbsPrShift );

                DMatrix& alpha = alphabeta->GetAlpha();
                DMatrix& beta  = alphabeta->GetBeta( );
            // ������ llr ���
                cllr->CalculateBitsMaxLLR( OutLLR, Alpha, Beta );

            // ������� �������� ������ �� NOFDM-��������
                SpectrsShift      += NL1; // ����� ������ �������� Spectrs
                InSymbsLogPrShift += N;   // InSymbsLogPr �
                BitNumberShift    += N*log2M; // OutLLR
        }

    // ��� ���������� � ������� mexFunction ������ (����� ��������
    // ����������) ���� ����������

        if( gamma ) 	delete gamma;
        if( alphabeta ) delete alphabeta;
        if( cllr )		delete cllr;
        

        delete [] SpectrsR;
        delete [] RefSigsR;
        if (isSigComplex){
            delete [] SpectrsI;
            delete [] RefSigsI;
        }

    return;
}

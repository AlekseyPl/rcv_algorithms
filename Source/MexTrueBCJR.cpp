#include "mex.h"
#include "math.h"

///////////////////////////////////////////////////////////////////////////
    void SetR2Zero(int NL1, int ML, double **R)
    // ������� ��������� ��������� R. Ÿ ���������� ��������� �����
    // �������� ���������� ���������� ���� ��������
    {
        int k, m;
        
        for (k = 0; k < NL1; k++)
            for (m = 0; m < ML; m++)
                *(R[k] + m) = 0;
    }

///////////////////////////////////////////////////////////////////////////
    void PreCalculateRWithAnyNt(int L, int N, int Nt, int ML, int ML1,
            int log2M, int SpectrsShift, double **Spectrs,
            double **RefSigs, double **R)
    // ������ ���������� ���������� ��� Nt > 1. ������� �� Nt = 1
    // ����������� � ������� ��������������� ����� for ��� �������
    // ����������
    {
        int k, m, n;
        int RefShift, SpectrsPos, MStep, MMax;
        double Buf;
        
        // ������ ��� ������ L-1 ��
            RefShift = 0;
            SpectrsPos = SpectrsShift;
            MStep = ML1;
            for (k = 0; k < L-1; k++) {
                for (m = 0; m < ML; m += MStep) {
                    for (n = 0; n < Nt; n++) {
                        Buf = *(RefSigs[RefShift + m] + n) -
                                *(Spectrs[SpectrsPos] + n);
                        *(R[k] + m) += Buf * Buf;
                    }
                }
                MStep >>= log2M;
                RefShift += ML;
                SpectrsPos++;
            }

        // ������ ��� ������� N-L+1 ��
            for (k = L-1; k < N; k++) {
                for (m = 0; m < ML; m++) {
                    for (n = 0; n < Nt; n++) {
                        Buf = *(RefSigs[RefShift + m] + n) -
                                *(Spectrs[SpectrsPos] + n);
                        *(R[k] + m) += Buf * Buf;
                    }
                }
                SpectrsPos++;
            }

        // ������ ��� ��������� L-1 ��
            RefShift += ML;
            MMax = ML1;
            for (k = N; k < N+L-1; k++) {
                for (m = 0; m < MMax; m++) {
                    for (n = 0; n < Nt; n++) {
                        Buf = *(RefSigs[RefShift + m] + n) -
                                *(Spectrs[SpectrsPos] + n);
                        *(R[k] + m) += Buf * Buf;
                    }
                }
                MMax >>= log2M;
                RefShift += ML;
                SpectrsPos++;
            }
        
        return;
    }

///////////////////////////////////////////////////////////////////////////
    void PreCalculateRWithOneNt(int L, int N, int ML, int ML1, int log2M, 
            int SpectrsShift, double **Spectrs, double **RefSigs,
            double **R)
    // ������ ���������� ���������� ��� Nt = 1. ������� �� Nt > 1
    // ����������� � ���������� ��������������� ����� for ��� �������
    // ����������
    {
        int k, m;
        int RefShift, SpectrsPos, MStep, MMax;
        double Buf;
        
        // ������ ��� ������ L-1 ��
            RefShift = 0;
            SpectrsPos = SpectrsShift;
            MStep = ML1;
            for (k = 0; k < L-1; k++) {
                for (m = 0; m < ML; m += MStep) {
                    Buf = *(RefSigs[RefShift + m]) -
                            *(Spectrs[SpectrsPos]);
                    *(R[k] + m) += Buf * Buf;
                }
                MStep >>= log2M;
                RefShift += ML;
                SpectrsPos++;
            }

        // ������ ��� ������� N-L+1 ��
            for (k = L-1; k < N; k++) {
                for (m = 0; m < ML; m++) {
                    Buf = *(RefSigs[RefShift + m]) -
                            *(Spectrs[SpectrsPos]);
                    *(R[k] + m) += Buf * Buf;
                }
                SpectrsPos++;
            }

        // ������ ��� ��������� L-1 ��
            RefShift += ML;
            MMax = ML1;
            for (k = N; k < N+L-1; k++) {
                for (m = 0; m < MMax; m++) {
                    Buf = *(RefSigs[RefShift + m]) -
                            *(Spectrs[SpectrsPos]);
                    *(R[k] + m) += Buf * Buf;
                }
                MMax >>= log2M;
                RefShift += ML;
                SpectrsPos++;
            }
        
        return;
    }

///////////////////////////////////////////////////////////////////////////
    void CalculateR(int L, int N, int Nt, int ML, int ML1, int log2M,
            int SpectrsShift, double Var, bool isSigComplex,
            double **SpectrsR, double **SpectrsI, double **RefSigsR,
            double **RefSigsI, double **R)
    // ������� ���������� ��������� �� ��������� ���������� ���������� �
    // ��������� ���������, ������ �� ������ �����
    {
        int k, m, MStep, MMax;
        
        // ������ �������� ���������� ����������
            if (Nt == 1) {
                PreCalculateRWithOneNt(    L, N,     ML, ML1, log2M,
                        SpectrsShift, SpectrsR, RefSigsR, R);
                if (isSigComplex)
                    PreCalculateRWithOneNt(L, N,     ML, ML1, log2M,
                            SpectrsShift, SpectrsI, RefSigsI, R);
            } else {
                PreCalculateRWithAnyNt(    L, N, Nt, ML, ML1, log2M,
                        SpectrsShift, SpectrsR, RefSigsR, R);
                if (isSigComplex)
                    PreCalculateRWithAnyNt(L, N, Nt, ML, ML1, log2M,
                            SpectrsShift, SpectrsI, RefSigsI, R);
            }

        // ���������� ��������� �� ���������� ���������� � ������ ��������
        // ��������� ����
            MStep = ML1;
            for (k = 0; k < L-1; k ++) {
                for (m = 0; m < ML; m += MStep)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // ��-������� � True Max MaxStar
                MStep >>= log2M;
            }

            for (k = L-1; k < N; k++) {
                for (m = 0; m < ML; m++)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // ��-������� � True Max MaxStar
            }

            MMax = ML1;
            for (k = N; k < N+L-1; k++) {
                for (m = 0; m < MMax; m++)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // ��-������� � True Max MaxStar
                MMax >>= log2M;
            }

        return;
    }

///////////////////////////////////////////////////////////////////////////
    void CalculateAlpha(int L, int N, int M, int ML1, int ML2, int log2M,
            int SymbsPrShift, double **Alpha, double **SymbsPr, double **R)
    {
        int n, t, SymbsPrState;
        int CurState, PrevState, RefSigNum, SymbNum;
        int CurStateStep, PrevStateStep, RefSigNumStep;
        int SymbNumCounter, SymbNumCounterMax;
        int PrevState0, MaxCurState;

        // ������ ������ ���������� � True Max MaxStar
        double Sum;
        
        
        
        
        // ����� ������ ���������� � True Max MaxStar
        
        CurStateStep  = ML2;
        PrevStateStep = ML1;
        RefSigNumStep = PrevStateStep;
        SymbNumCounterMax = 1;
        SymbsPrState = SymbsPrShift;
        for (t = 1; t < L; t++){
            PrevState      = 0;
            RefSigNum      = 0;
            SymbNum        = 0;
            SymbNumCounter = 0;
            Sum            = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState += CurStateStep){
                // ������ ������ ���������� � True Max MaxStar
                *(Alpha[t] + CurState) = 
                        *(Alpha[t-1] + PrevState) * 
                        *(R[t-1] + RefSigNum) *
                        *(SymbsPr[SymbNum]+ SymbsPrState);
                Sum += *(Alpha[t] + CurState);
                
                // ����� ������ ���������� � True Max MaxStar
                PrevState += PrevStateStep;
                if (PrevState >= ML1)
                    PrevState -= ML1;
                RefSigNum += RefSigNumStep;
                SymbNumCounter++;
                if (SymbNumCounter == SymbNumCounterMax){
                    SymbNum++;
                    SymbNumCounter = 0;
                }
            }
            
            for (CurState = 0; CurState < ML1; CurState += CurStateStep)
                *(Alpha[t] + CurState) /= Sum; // ��-������� � True Max MaxStar

            PrevStateStep = CurStateStep;
            RefSigNumStep = PrevStateStep;
            CurStateStep >>= log2M;
            SymbNumCounterMax <<= log2M;
            SymbsPrState++;
        }

        for (t = L; t <= N; t++){
            PrevState0     = 0;
            RefSigNum      = 0;
            SymbNum        = 0;
            SymbNumCounter = 0;
            Sum            = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState++){
                *(Alpha[t] + CurState) = 0;  // ��-������� � True Max MaxStar
                for (n = 0; n < M; n++){
                    // ������ ������ ���������� � True Max MaxStar
                    *(Alpha[t] + CurState) += 
                        *(Alpha[t-1] + PrevState0 + n) * 
                        *(R[t-1] + RefSigNum) *
                        *(SymbsPr[SymbNum]+ SymbsPrState);


                    // ����� ������ ���������� � True Max MaxStar
                    RefSigNum++;
                }
                // ������ ������ ���������� � True Max MaxStar
                Sum += *(Alpha[t] + CurState);
                
                // ����� ������ ���������� � True Max MaxStar
                PrevState0 += M;
                if (PrevState0 >= ML1)
                    PrevState0 -= ML1;
                SymbNumCounter++;
                if (SymbNumCounter == ML2){
                    SymbNum++;
                    SymbNumCounter = 0;
                }
            }

            for (CurState = 0; CurState < ML1; CurState++)
                *(Alpha[t] + CurState) /= Sum; // ��-������� � True Max MaxStar
            
            SymbsPrState++;
        }

        MaxCurState = ML2;
        for (t = N+1; t < N+L; t++){
            PrevState0 = 0;
            RefSigNum  = 0;
            Sum        = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < MaxCurState; CurState++){
                *(Alpha[t] + CurState) = 0; // ��-������� � True Max MaxStar
                for (n = 0; n < M; n++){
                    // ������ ������ ���������� � True Max MaxStar
                    *(Alpha[t] + CurState) += 
                        *(Alpha[t-1] + PrevState0 + n) * 
                        *(R[t-1] + RefSigNum);
                    
                    
                    // ����� ������ ���������� � True Max MaxStar
                    RefSigNum++;
                }
                
                Sum += *(Alpha[t] + CurState);
                PrevState0 += M;
            }

            for (CurState = 0; CurState < MaxCurState; CurState++)
                *(Alpha[t] + CurState) /= Sum; // ��-������� � True Max MaxStar

            MaxCurState >>= log2M;
        }
        
        return;
    }

///////////////////////////////////////////////////////////////////////////
    void CalculateBeta(int L, int N, int M, int ML1, int ML2, int log2M,
            int SymbsPrShift, double **Beta, double **SymbsPr, double **R)
    {
        int n, t, SymbsPrState;
        int CurState, NextState, NextState0, RefSigNum, RefSigNum0;
        int MaxCurState, NextState0Step, RefSigNum0Step;
        int NextStateCount, NextState0Count;

        // ������ ������ ���������� � True Max MaxStar
        double Sum;
        
        
        
        
        // ����� ������ ���������� � True Max MaxStar
        
        MaxCurState = M;
        for (t = N+L-2; t >= N; t--) {
            NextState      = 0;
            NextStateCount = 0;
            RefSigNum      = 0;
            Sum            = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < MaxCurState; CurState++) {
                // ������ ������ ���������� � True Max MaxStar
                *(Beta[t] + CurState) = 
                        *(Beta[t+1] + NextState) *
                        *(R[t] + RefSigNum);
                Sum += *(Beta[t] + CurState);
                
                // ����� ������ ���������� � True Max MaxStar
                NextStateCount++;
                if (NextStateCount == M) {
                    NextStateCount = 0;
                    NextState++;
                }
                RefSigNum++;
            }

            for (CurState = 0; CurState < MaxCurState; CurState++)
                *(Beta[t] + CurState) /= Sum; // ��-������� � True Max MaxStar

            MaxCurState <<= log2M;
        }
    
        NextState0Count = 0;
        SymbsPrState = SymbsPrShift + N-1;
        for (t = N-1; t >= L-1; t--) {
            NextState0 = 0;
            RefSigNum0 = 0;
            Sum        = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState++) {
                NextState = NextState0;
                RefSigNum = RefSigNum0;
                *(Beta[t] + CurState) = 0; // ��-������� � True Max MaxStar
                for (n = 0; n < M; n++) {
                    // ������ ������ ���������� � True Max MaxStar
                    *(Beta[t] + CurState) += 
                            *(Beta[t+1] + NextState) * 
                            *(R[t] + RefSigNum) *
                            *(SymbsPr[n]+ SymbsPrState);
                    
                    
                    // ����� ������ ���������� � True Max MaxStar
                    NextState += ML2;
                    RefSigNum += ML1;
                }
                // ������ ������ ���������� � True Max MaxStar
                Sum += *(Beta[t] + CurState);
                
                // ����� ������ ���������� � True Max MaxStar
                
                NextState0Count++;
                if (NextState0Count == M) {
                    NextState0Count = 0;
                    NextState0++;
                }
                
                RefSigNum0++;
            }

            for (CurState = 0; CurState < ML1; CurState++)
                *(Beta[t] + CurState) /= Sum; // ��-������� � True Max MaxStar
            
            SymbsPrState--;
        }

        NextState0Step = 1;
        RefSigNum0Step = M;
        for (t = L-2; t >= 0; t--){
            NextState0 = 0;
            RefSigNum0 = 0;
            Sum        = 0; // ��-������� � True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState += RefSigNum0Step) {
                NextState = NextState0;
                RefSigNum = RefSigNum0;
                *(Beta[t] + CurState) = 0; // ��-������� � True Max MaxStar
                for (n = 0; n < M; n++) {
                    // ������ ������ ���������� � True Max MaxStar
                    *(Beta[t] + CurState) += 
                            *(Beta[t+1] + NextState) * 
                            *(R[t] + RefSigNum) *
                            *(SymbsPr[n] + SymbsPrState);
                    
                    
                    // ����� ������ ���������� � True Max MaxStar
                    NextState += ML2;
                    RefSigNum += ML1;
                }
                // ������ ������ ���������� � True Max MaxStar
                Sum += *(Beta[t] + CurState);
                
                // ����� ������ ���������� � True Max MaxStar
                RefSigNum0 += RefSigNum0Step;
                NextState0 += NextState0Step;
            }

            for (CurState = 0; CurState < ML1; CurState += RefSigNum0Step)
                *(Beta[t] + CurState) /= Sum; // ��-������� � True Max MaxStar
            
            NextState0Step = RefSigNum0Step;
            RefSigNum0Step <<= log2M;
            SymbsPrState--;
        }

        return;
    }

///////////////////////////////////////////////////////////////////////////
//     void CalculateSymbsPr(int N, int M, int ML2, int SymbsPrShift,
//             double **SymbsPr, double **Alpha, double **Beta)
//     {
//         int k, m, t;
//         int StateShift, SymbsPrState;
//         
//         StateShift = 0;
//         for (m = 0; m < M; m++) {
//             SymbsPrState = SymbsPrShift;
//             for (t = 1; t < N+1; t++) {
//                 *(SymbsPr[m] + SymbsPrState) = 0;
//                 for (k = StateShift; k < StateShift + ML2; k++) {
//                     *(SymbsPr[m] + SymbsPrState) +=
//                             *(Alpha[t] + k) * *(Beta[t] + k);
//                 }
//                 SymbsPrState++;
//             }
//             StateShift += ML2;
//         }
// 
//         return;
//     }

///////////////////////////////////////////////////////////////////////////
    void CalculateBitsLLR(int L, int N, int M, int ML, int ML1, int log2M,
            int BitNumberShift, double *OutLLR, double **Alpha,
            double **Beta, double **Lambda)
    {
        int k, m, t;
        int BitNumber, NumStates, StateStep, StateShift;
        // ������ ������ ���������� � True Max MaxStar
        double PrZer, PrOne;
        
        
        
        // ����� ������ ���������� � True Max MaxStar
        
        for (t = 1; t < N+1; t++)
            for (k = 0; k < ML1; k++)
                *(Lambda[t] + k) = *(Alpha[t] + k) * *(Beta[t] + k); // ��-������� � True Max MaxStar
        
        BitNumber = BitNumberShift;
        StateStep = ML1 >> log2M;
        for (t = 1; t < L-1; t++) {
            NumStates = (ML1 >> 1);
            for (m = 0; m < log2M; m++) {
                PrZer = 0; // ��-������� � True Max MaxStar
                PrOne = 0; // ��-������� � True Max MaxStar
                StateShift = 0;
                while (StateShift < ML1){
                    for (k = StateShift; k < StateShift + NumStates; k += StateStep)
                        PrZer += *(Lambda[t] + k); // ��-������� � True Max MaxStar
                    
                    StateShift += NumStates;
                    for (k = StateShift; k < StateShift + NumStates; k += StateStep)
                        PrOne += *(Lambda[t] + k); // ��-������� � True Max MaxStar
                    
                    StateShift += NumStates;
                }
                *(OutLLR + BitNumber) = log(PrZer / PrOne); // ��-������� � True Max MaxStar
                BitNumber++;
                NumStates >>= 1;
            }
            StateStep >>= log2M;
        }

        for (t = L-1; t < N+1; t++) {
            NumStates = (ML1 >> 1);
            for (m = 0; m < log2M; m++) {
                PrZer = 0; // ��-������� � True Max MaxStar
                PrOne = 0; // ��-������� � True Max MaxStar
                StateShift = 0;
                while (StateShift < ML1){
                    for (k = StateShift; k < StateShift + NumStates; k++)
                        PrZer += *(Lambda[t] + k); // ��-������� � True Max MaxStar
                    
                    StateShift += NumStates;
                    for (k = StateShift; k < StateShift + NumStates; k++)
                        PrOne += *(Lambda[t] + k); // ��-������� � True Max MaxStar
                    
                    StateShift += NumStates;
                }
                *(OutLLR + BitNumber) = log(PrZer / PrOne); // ��-������� � True Max MaxStar
                BitNumber++;
                NumStates >>= 1;
            }
        }
        
        return;
    }
    
///////////////////////////////////////////////////////////////////////////
    void LLR2SymbsPr(int M, int log2M, int NNSymb, double *LLR,
            double **SymbsPr)
    {
        int k, n, m, Mask;
        int **Map;
        double **BitsPr;

        // ����������� llr � ����������� ��������� 0 � 1 �
        // ��������� �� � ������� BitsPr �������� N*NSymb x 2*log2M
            // ��������� �� ������ BitsPr
                BitsPr = new double*[log2M+log2M];
                for (k = 0; k < log2M+log2M; k++)
                    BitsPr[k] = new double[NNSymb];
        
            m = 0;
            for (k = 0; k < NNSymb; k++)
                for (n = 0; n < log2M; n++){
                // ������ ������ ���������� � True Max MaxStar
                    // ����������� ��������� 1
                        *(BitsPr[n + log2M] + k) = 
                                1 / (1 + exp(*(LLR + m)));
                        
                        
                    // ����������� ��������� 0
                        *(BitsPr[n] + k) = 1 -
                                *(BitsPr[n + log2M] + k);
                    // �������������� ������ � ������� llr
                        m++;
                // ����� ������ ���������� � True Max MaxStar
                }

        // �������� ����� ������ ��������� BitsPr ��� ����������
        // ������������ ��������. ������ ����� M x log2M
            // ��������� �� ������ Map
                Map = new int*[log2M];
                for (k = 0; k < log2M; k++)
                    Map[k] = new int[M];
            
            for (k = 0; k < M; k++){
                Mask = 1 << log2M;
                for (n = 0; n < log2M; n++){
                    Mask >>= 1;
                    if (n && Mask == 0)
                        *(Map[n] + k) = n;
                    else
                        *(Map[n] + k) = n + log2M;
                }
            }

        // �������� ����������� �������� � ������� �� � ������ SymbsPr
        // �������� N*NSymb x M
            for (k = 0; k < NNSymb; k++)
                for (n = 0; n < M; n++){
                // ������ ������ ���������� � True Max MaxStar
                    *(SymbsPr[n] + k) = 1;
                    for (m = 0; m < log2M; m++)
                        *(SymbsPr[n] + k) *=
                                *(BitsPr[*(Map[m] + n)] + k);
                // ����� ������ ���������� � True Max MaxStar
                }

        // ������������ ���������� ������
            for (k = 0; k < log2M+log2M; k++)
                delete [] BitsPr[k];
            delete [] BitsPr;

            for (k = 0; k < log2M; k++)
                delete [] Map[k];
            delete [] Map;

        return;
    }

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // ��������! � ������� ��� ������� ������ �� ������, ������� ���
    // "���������" ������� ������ ����� ���� ��������� ���� � �������
    // ������� ����������!
    //
    // ������ ������ mex-������� �� MATLAB:
    //
    //   OutLLR = MexTrueBCJR(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
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
            // ��������� �� ������ SymbsPr
                double *OutLLR;

        // ���������� ����������
            // ����� ������������ ��������: log2M = log2(M), ML = M^L;
            // ML1 = M^(L-1); ML2 = M^(L-2); NL = N+L; NL1 = N+L-1,
            // NNSymb = N*NSymb
                int log2M, ML, ML1, ML2, NL, NL1, NNSymb;
            // ��������
                int k, SpectrsShift, InSymbsPrShift, BitNumberShift;
            // ����, ����������� �� �� ����������� ������ ��� ���
                bool isSigComplex;
            // R, Alpha, Beta
                double **R, **Alpha, **Beta, **Lambda, **InSymbsPr;

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
        // ��������� �� ������ Spectr
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

        // ��������� �� ������ InLLR
            InLLR = mxGetPr(prhs[8]);

    // ��������� ������ ��� �������� ���������� - ������ OutLLR. ���
	// ��������� ������ ���� ������ ������ ��������, ����� MATLAB ����
    // �������� ���������.
        plhs[0] = mxCreateDoubleMatrix(NNSymb*log2M, 1, mxREAL);
        OutLLR = mxGetPr(plhs[0]);
            
    // ��������� ������ ��� ���������� R, Alpha, Beta, Lambda, InSymbsPr
    // ���� �� ����� ��������� ���� ��������, �� ��� ��� ����������
    // ����������, ����� ��� ����������� ��� �������� ����������
    // ��� ��������� ������ ����� double ���������� ���������� ������ �����
    // ���� �����, ��� ��������� ����� mxCreateDoubleMatrix ������������
    // ���������� ������! �������� �������� ���, ��� R ���� �������� ������
    // ��� ������� �������; Alpha, Beta - ����� ���������������� ������
    // ���� �������� � ������ ���� ���; Lambda - �� ����� ����������������;
    // InSymbsPr - �������� ���������������� ������ LLR2SymbsPr
        if (nlhs == 1) {
            // ��������� �� ������ R (������ ML x NL1)
                R = new double*[NL1];
                for (k = 0; k < NL1; k++)
                    R[k] = new double[ML];

            // ������������� ���������� �� ������� Alpha, Beta � Lambda
            // (������ ML1 x NL), ������� ������������� �������� Alpha,
            // Beta
                Alpha = new double*[NL];
                for (k = 0; k < NL; k++)
                    Alpha[k] = new double[ML1];
                *(Alpha[0] + 0) = 1; // ��-������� � True Max MaxStar

                Beta = new double*[NL];
                for (k = 0; k < NL; k++)
                    Beta[k] = new double[ML1];
                *(Beta[NL1] + 0) = 1; // ��-������� � True Max MaxStar

                Lambda = new double*[NL];
                for (k = 0; k < NL; k++)
                    Lambda[k] = new double[ML1];

            // ��������� �� ������ InSymbsPr (������ N*NSymb x M)
                InSymbsPr = new double*[M];
                for (k = 0; k < M; k++)
                    InSymbsPr[k] = new double[NNSymb];
        } else {
            plhs[1] = mxCreateDoubleMatrix(ML, NL1, mxREAL);
            R = new double*[NL1];
            R[0] = mxGetPr(plhs[1]);
            for (k = 1; k < NL1; k++)
                R[k] = R[k-1] + ML;

            plhs[2] = mxCreateDoubleMatrix(ML1, NL, mxREAL);
            Alpha = new double*[NL];
            Alpha[0] = mxGetPr(plhs[2]);
            for (k = 1; k < NL; k++)
                Alpha[k] = Alpha[k-1] + ML1;
            *(Alpha[0] + 0) = 1; // ��-������� � True Max MaxStar

            plhs[3] = mxCreateDoubleMatrix(ML1, NL, mxREAL);
            Beta = new double*[NL];
            Beta[0] = mxGetPr(plhs[3]);
            for (k = 1; k < NL; k++)
                Beta[k] = Beta[k-1] + ML1;
            *(Beta[NL1] + 0) = 1; // ��-������� � True Max MaxStar

            plhs[4] = mxCreateDoubleMatrix(ML1, NL, mxREAL);
            Lambda = new double*[NL];
            Lambda[0] = mxGetPr(plhs[4]);
            for (k = 1; k < NL; k++)
                Lambda[k] = Lambda[k-1] + ML1;

            plhs[5] = mxCreateDoubleMatrix(NNSymb, M, mxREAL);
            InSymbsPr = new double*[M];
            InSymbsPr[0] = mxGetPr(plhs[5]);
            for (k = 1; k < M; k++)
                InSymbsPr[k] = InSymbsPr[k-1] + NNSymb;
        }

    // ������ ���� ��������� - �������������� ������� llr � �����������
    // �������� ��� ����� �����
        LLR2SymbsPr(M, log2M, NNSymb, InLLR, InSymbsPr);

    // ���� �� ���������� NOFDM-��������
        SpectrsShift   = 0; // ����� ������ �������� Spectrs,
        InSymbsPrShift = 0; // InSymbsPr �
        BitNumberShift = 0; // OutLLR
        for (k = 0; k < NSymb; k++) {
            // ������ R - �������� ������������ �������� ������� �������� �
            // �������� �������
                SetR2Zero(NL1, ML, R);
                CalculateR(L, N, Nt, ML, ML1, log2M, SpectrsShift, Var,
                        isSigComplex, SpectrsR, SpectrsI, RefSigsR,
                        RefSigsI, R);

            // ������ ������� Alpha
                CalculateAlpha(L, N, M, ML1, ML2, log2M, InSymbsPrShift,
                        Alpha, InSymbsPr, R);

            // ������ ������� Beta
                CalculateBeta (L, N, M, ML1, ML2, log2M, InSymbsPrShift,
                        Beta, InSymbsPr, R);

            // ������ llr ���
                CalculateBitsLLR(L, N, M, ML, ML1, log2M, BitNumberShift,
                        OutLLR, Alpha, Beta, Lambda);

            // ������� �������� ������ �� NOFDM-��������
                SpectrsShift   += NL1; // ����� ������ �������� Spectrs,
                InSymbsPrShift += N;   // InSymbsPr �
                BitNumberShift += N*log2M; // OutLLR
        }

    // ��� ���������� � ������� mexFunction ������ (����� ��������
    // ����������) ���� ����������
        if (nlhs == 1) {
            for (k = 0; k < NL1; k++)
                delete [] R[k];
            for (k = 0; k < NL; k++)
                delete [] Alpha[k];
            for (k = 0; k < NL; k++)
                delete [] Beta[k];
            for (k = 0; k < NL; k++)
                delete [] Lambda[k];
            for (k = 0; k < M; k++)
                delete [] InSymbsPr[k];
        }
        
        delete [] R;
        delete [] Alpha;
        delete [] Beta;
        delete [] Lambda;
        delete [] InSymbsPr;

        delete [] SpectrsR;
        delete [] RefSigsR;
        if (isSigComplex){
            delete [] SpectrsI;
            delete [] RefSigsI;
        }
        
    return;
}
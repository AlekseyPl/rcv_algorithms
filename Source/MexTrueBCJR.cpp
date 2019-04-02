#include "mex.h"
#include "math.h"

///////////////////////////////////////////////////////////////////////////
    void SetR2Zero(int NL1, int ML, double **R)
    // Функция обнуления элементов R. Её необходимо выполнять перед
    // расчётом евклидовых расстояний всех символов
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
    // Расчёт евклидовых расстояний для Nt > 1. Отличие от Nt = 1
    // заключается в наличие дополнительного цикла for для каждого
    // расстояния
    {
        int k, m, n;
        int RefShift, SpectrsPos, MStep, MMax;
        double Buf;
        
        // Расчёт для первых L-1 ТИ
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

        // Расчёт для средних N-L+1 ТИ
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

        // Расчёт для последних L-1 ТИ
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
    // Расчёт евклидовых расстояний для Nt = 1. Отличие от Nt > 1
    // заключается в отсутствии дополнительного цикла for для каждого
    // расстояния
    {
        int k, m;
        int RefShift, SpectrsPos, MStep, MMax;
        double Buf;
        
        // Расчёт для первых L-1 ТИ
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

        // Расчёт для средних N-L+1 ТИ
            for (k = L-1; k < N; k++) {
                for (m = 0; m < ML; m++) {
                    Buf = *(RefSigs[RefShift + m]) -
                            *(Spectrs[SpectrsPos]);
                    *(R[k] + m) += Buf * Buf;
                }
                SpectrsPos++;
            }

        // Расчёт для последних L-1 ТИ
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
    // Функция вычисления экспонент от отношений эвклидовых расстояний к
    // удвоенной дисперсии, взятых со знаком минус
    {
        int k, m, MStep, MMax;
        
        // Расчёт значений эвклидовых расстояний
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

        // Вычисление экспонент от эвклидовых расстояний с учётом значения
        // дисперсии шума
            MStep = ML1;
            for (k = 0; k < L-1; k ++) {
                for (m = 0; m < ML; m += MStep)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // По-разному в True Max MaxStar
                MStep >>= log2M;
            }

            for (k = L-1; k < N; k++) {
                for (m = 0; m < ML; m++)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // По-разному в True Max MaxStar
            }

            MMax = ML1;
            for (k = N; k < N+L-1; k++) {
                for (m = 0; m < MMax; m++)
                    *(R[k] + m) = exp(*(R[k] + m) * Var); // По-разному в True Max MaxStar
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

        // Начало разной реализации в True Max MaxStar
        double Sum;
        
        
        
        
        // Конец разной реализации в True Max MaxStar
        
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
            Sum            = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState += CurStateStep){
                // Начало разной реализации в True Max MaxStar
                *(Alpha[t] + CurState) = 
                        *(Alpha[t-1] + PrevState) * 
                        *(R[t-1] + RefSigNum) *
                        *(SymbsPr[SymbNum]+ SymbsPrState);
                Sum += *(Alpha[t] + CurState);
                
                // Конец разной реализации в True Max MaxStar
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
                *(Alpha[t] + CurState) /= Sum; // По-разному в True Max MaxStar

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
            Sum            = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState++){
                *(Alpha[t] + CurState) = 0;  // По-разному в True Max MaxStar
                for (n = 0; n < M; n++){
                    // Начало разной реализации в True Max MaxStar
                    *(Alpha[t] + CurState) += 
                        *(Alpha[t-1] + PrevState0 + n) * 
                        *(R[t-1] + RefSigNum) *
                        *(SymbsPr[SymbNum]+ SymbsPrState);


                    // Конец разной реализации в True Max MaxStar
                    RefSigNum++;
                }
                // Начало разной реализации в True Max MaxStar
                Sum += *(Alpha[t] + CurState);
                
                // Конец разной реализации в True Max MaxStar
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
                *(Alpha[t] + CurState) /= Sum; // По-разному в True Max MaxStar
            
            SymbsPrState++;
        }

        MaxCurState = ML2;
        for (t = N+1; t < N+L; t++){
            PrevState0 = 0;
            RefSigNum  = 0;
            Sum        = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < MaxCurState; CurState++){
                *(Alpha[t] + CurState) = 0; // По-разному в True Max MaxStar
                for (n = 0; n < M; n++){
                    // Начало разной реализации в True Max MaxStar
                    *(Alpha[t] + CurState) += 
                        *(Alpha[t-1] + PrevState0 + n) * 
                        *(R[t-1] + RefSigNum);
                    
                    
                    // Конец разной реализации в True Max MaxStar
                    RefSigNum++;
                }
                
                Sum += *(Alpha[t] + CurState);
                PrevState0 += M;
            }

            for (CurState = 0; CurState < MaxCurState; CurState++)
                *(Alpha[t] + CurState) /= Sum; // По-разному в True Max MaxStar

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

        // Начало разной реализации в True Max MaxStar
        double Sum;
        
        
        
        
        // Конец разной реализации в True Max MaxStar
        
        MaxCurState = M;
        for (t = N+L-2; t >= N; t--) {
            NextState      = 0;
            NextStateCount = 0;
            RefSigNum      = 0;
            Sum            = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < MaxCurState; CurState++) {
                // Начало разной реализации в True Max MaxStar
                *(Beta[t] + CurState) = 
                        *(Beta[t+1] + NextState) *
                        *(R[t] + RefSigNum);
                Sum += *(Beta[t] + CurState);
                
                // Конец разной реализации в True Max MaxStar
                NextStateCount++;
                if (NextStateCount == M) {
                    NextStateCount = 0;
                    NextState++;
                }
                RefSigNum++;
            }

            for (CurState = 0; CurState < MaxCurState; CurState++)
                *(Beta[t] + CurState) /= Sum; // По-разному в True Max MaxStar

            MaxCurState <<= log2M;
        }
    
        NextState0Count = 0;
        SymbsPrState = SymbsPrShift + N-1;
        for (t = N-1; t >= L-1; t--) {
            NextState0 = 0;
            RefSigNum0 = 0;
            Sum        = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState++) {
                NextState = NextState0;
                RefSigNum = RefSigNum0;
                *(Beta[t] + CurState) = 0; // По-разному в True Max MaxStar
                for (n = 0; n < M; n++) {
                    // Начало разной реализации в True Max MaxStar
                    *(Beta[t] + CurState) += 
                            *(Beta[t+1] + NextState) * 
                            *(R[t] + RefSigNum) *
                            *(SymbsPr[n]+ SymbsPrState);
                    
                    
                    // Конец разной реализации в True Max MaxStar
                    NextState += ML2;
                    RefSigNum += ML1;
                }
                // Начало разной реализации в True Max MaxStar
                Sum += *(Beta[t] + CurState);
                
                // Конец разной реализации в True Max MaxStar
                
                NextState0Count++;
                if (NextState0Count == M) {
                    NextState0Count = 0;
                    NextState0++;
                }
                
                RefSigNum0++;
            }

            for (CurState = 0; CurState < ML1; CurState++)
                *(Beta[t] + CurState) /= Sum; // По-разному в True Max MaxStar
            
            SymbsPrState--;
        }

        NextState0Step = 1;
        RefSigNum0Step = M;
        for (t = L-2; t >= 0; t--){
            NextState0 = 0;
            RefSigNum0 = 0;
            Sum        = 0; // По-разному в True Max MaxStar
            for (CurState = 0; CurState < ML1; CurState += RefSigNum0Step) {
                NextState = NextState0;
                RefSigNum = RefSigNum0;
                *(Beta[t] + CurState) = 0; // По-разному в True Max MaxStar
                for (n = 0; n < M; n++) {
                    // Начало разной реализации в True Max MaxStar
                    *(Beta[t] + CurState) += 
                            *(Beta[t+1] + NextState) * 
                            *(R[t] + RefSigNum) *
                            *(SymbsPr[n] + SymbsPrState);
                    
                    
                    // Конец разной реализации в True Max MaxStar
                    NextState += ML2;
                    RefSigNum += ML1;
                }
                // Начало разной реализации в True Max MaxStar
                Sum += *(Beta[t] + CurState);
                
                // Конец разной реализации в True Max MaxStar
                RefSigNum0 += RefSigNum0Step;
                NextState0 += NextState0Step;
            }

            for (CurState = 0; CurState < ML1; CurState += RefSigNum0Step)
                *(Beta[t] + CurState) /= Sum; // По-разному в True Max MaxStar
            
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
        // Начало разной реализации в True Max MaxStar
        double PrZer, PrOne;
        
        
        
        // Конец разной реализации в True Max MaxStar
        
        for (t = 1; t < N+1; t++)
            for (k = 0; k < ML1; k++)
                *(Lambda[t] + k) = *(Alpha[t] + k) * *(Beta[t] + k); // По-разному в True Max MaxStar
        
        BitNumber = BitNumberShift;
        StateStep = ML1 >> log2M;
        for (t = 1; t < L-1; t++) {
            NumStates = (ML1 >> 1);
            for (m = 0; m < log2M; m++) {
                PrZer = 0; // По-разному в True Max MaxStar
                PrOne = 0; // По-разному в True Max MaxStar
                StateShift = 0;
                while (StateShift < ML1){
                    for (k = StateShift; k < StateShift + NumStates; k += StateStep)
                        PrZer += *(Lambda[t] + k); // По-разному в True Max MaxStar
                    
                    StateShift += NumStates;
                    for (k = StateShift; k < StateShift + NumStates; k += StateStep)
                        PrOne += *(Lambda[t] + k); // По-разному в True Max MaxStar
                    
                    StateShift += NumStates;
                }
                *(OutLLR + BitNumber) = log(PrZer / PrOne); // По-разному в True Max MaxStar
                BitNumber++;
                NumStates >>= 1;
            }
            StateStep >>= log2M;
        }

        for (t = L-1; t < N+1; t++) {
            NumStates = (ML1 >> 1);
            for (m = 0; m < log2M; m++) {
                PrZer = 0; // По-разному в True Max MaxStar
                PrOne = 0; // По-разному в True Max MaxStar
                StateShift = 0;
                while (StateShift < ML1){
                    for (k = StateShift; k < StateShift + NumStates; k++)
                        PrZer += *(Lambda[t] + k); // По-разному в True Max MaxStar
                    
                    StateShift += NumStates;
                    for (k = StateShift; k < StateShift + NumStates; k++)
                        PrOne += *(Lambda[t] + k); // По-разному в True Max MaxStar
                    
                    StateShift += NumStates;
                }
                *(OutLLR + BitNumber) = log(PrZer / PrOne); // По-разному в True Max MaxStar
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

        // Преобразуем llr в вероятности появления 0 и 1 и
        // расставим их в массиве BitsPr размером N*NSymb x 2*log2M
            // Указатель на массив BitsPr
                BitsPr = new double*[log2M+log2M];
                for (k = 0; k < log2M+log2M; k++)
                    BitsPr[k] = new double[NNSymb];
        
            m = 0;
            for (k = 0; k < NNSymb; k++)
                for (n = 0; n < log2M; n++){
                // Начало разной реализации в True Max MaxStar
                    // Вероятности появления 1
                        *(BitsPr[n + log2M] + k) = 
                                1 / (1 + exp(*(LLR + m)));
                        
                        
                    // Вероятности появления 0
                        *(BitsPr[n] + k) = 1 -
                                *(BitsPr[n + log2M] + k);
                    // Инкрементируем индекс в массиве llr
                        m++;
                // Конец разной реализации в True Max MaxStar
                }

        // Составим карту выбора элементов BitsPr для вычисления
        // вероятностей символов. Размер карты M x log2M
            // Указатель на массив Map
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

        // Вычислим вероятности символов и запишем их в массив SymbsPr
        // размером N*NSymb x M
            for (k = 0; k < NNSymb; k++)
                for (n = 0; n < M; n++){
                // Начало разной реализации в True Max MaxStar
                    *(SymbsPr[n] + k) = 1;
                    for (m = 0; m < log2M; m++)
                        *(SymbsPr[n] + k) *=
                                *(BitsPr[*(Map[m] + n)] + k);
                // Конец разной реализации в True Max MaxStar
                }

        // Освобождение выделенной памяти
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
    // Внимание! В функции нет никакой защиты от дурака, поэтому при
    // "вылетании" функции прежде всего надо проверять типы и размеры
    // входных переменных!
    //
    // Пример вызова mex-функции из MATLAB:
    //
    //   OutLLR = MexTrueBCJR(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
    //          InLLR, TBLen);
    //
    //   Входные переменные:
    //      L     - кол-во компонент, натуральное число;
    //      M     - размер созвездия, натуральное число;
    //      N     - кол-во модуляционных символов в одном NOFDM-символе,
    //              натуральное число;
    //      Nt    - кол-во отсчётов сигнала на один ТИ, натуральное число;
    //      NSymb - кол-во NOFDM-символов в одном кадре, натуральное число;
    //      Var   - дисперсия АБГШ, вещественное число;
    //      Spectrs - (Nt*(N+L-1) x NSymb) массив вещественных /
    //                комплексных спектров - здесь он будет трактоваться
    //                как массив (Nt x (N+L-1)*NSymb);
    //      RefSigs - (Nt х (M^L)*(2*L-1)) массив вещественных /
    //                комплексных опорных сигналов;
    //      InLLR   - (N*NSymb*log2(M) х 1) массив вещественных значений
    //                априорных llr бит;
    //      TBLen   - возможные значения от 1 до (N+L-1) - НЕ ИСПОЛЬЗУЕТСЯ!
    //
    //   Выходные переменные:
    //      OutLLR  - (N*NSymb*log2(M) х 1) массив вещественных значений
    //                апостериорных llr бит.
    //
    // Внимание к матричным переменным! При передаче двумерного массива,
    // т.е. матрицы, из Matlab в C++ данные записываются в одномерный
    // массив, причём порядок считывания из матрицы - по столбцам!

    // Объявления переменных
        // Входные переменные
            // целочисленные переменные-скаляры
                int L, M, N, Nt, NSymb;
                int TBLen;
            // вещественные переменные-скаляры
                double Var;
            // Указатели на вещественную (и мнимую) часть(и)
            // переменных-массивов
                double **SpectrsR, **SpectrsI, **RefSigsR, **RefSigsI,
                        *InLLR;

        // Выходные переменные
            // Указатель на массив SymbsPr
                double *OutLLR;

        // Внутренние переменные
            // Часто используемые значения: log2M = log2(M), ML = M^L;
            // ML1 = M^(L-1); ML2 = M^(L-2); NL = N+L; NL1 = N+L-1,
            // NNSymb = N*NSymb
                int log2M, ML, ML1, ML2, NL, NL1, NNSymb;
            // Счётчики
                int k, SpectrsShift, InSymbsPrShift, BitNumberShift;
            // Флаг, указывающий на то комплексный сигнал или нет
                bool isSigComplex;
            // R, Alpha, Beta
                double **R, **Alpha, **Beta, **Lambda, **InSymbsPr;

    // Получение значений входных переменных
        // Количество компонент
            L  = (int)(mxGetScalar(prhs[0]));
        // Размер созвездия
            M  = (int)(mxGetScalar(prhs[1]));
        // Кол-во модуляционных символов в одном NOFDM-символе 
            N  = (int)(mxGetScalar(prhs[2]));
        // Кол-во отсчётов сигнала на один ТИ
            Nt = (int)(mxGetScalar(prhs[3]));
        // Кол-во NOFDM-символов в одном кадре
            NSymb = (int)(mxGetScalar(prhs[4]));
        // Дисперсия шума
            Var = (double)(mxGetScalar(prhs[5]));
        // TBLen - НЕ ИСПОЛЬЗУЕТСЯ!
        //    TBLen = (int)(mxGetScalar(prhs[9]));
            
    // Расчёт часто используемых значений - для упрощения кода    
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
        // -1/(2*Var) - множитель показателя экспоненты
            Var = -1/(2*Var);

        // Флаг, указывающий на то комплексный сигнал или нет
            isSigComplex = mxIsComplex(prhs[6]);
            
    // Продолжение получения значений входных переменных
        // Указатели на массив Spectr
            // Вещественная часть
                SpectrsR = new double*[NL1*NSymb];
                SpectrsR[0] = mxGetPr(prhs[6]);
                for (k = 1; k < NL1*NSymb; k++)
                    SpectrsR[k] = SpectrsR[k-1] + Nt;
            // Мнимая часть
                if (isSigComplex) {
                    SpectrsI = new double*[NL1*NSymb];
                    SpectrsI[0] = mxGetPi(prhs[6]);
                    for (k = 1; k < NL1*NSymb; k++)
                        SpectrsI[k] = SpectrsI[k-1] + Nt;
                }
            
        // Указатели на массив RefSigs
            // Вещественная часть
                RefSigsR = new double*[ML*(2*L-1)];
                RefSigsR[0] = mxGetPr(prhs[7]);
                for (k = 1; k < ML*(2*L-1); k++)
                    RefSigsR[k] = RefSigsR[k-1] + Nt;
            // Мнимая часть
                if (isSigComplex) {                
                    RefSigsI = new double*[ML*(2*L-1)];
                    RefSigsI[0] = mxGetPi(prhs[7]);
                    for (k = 1; k < ML*(2*L-1); k++)
                        RefSigsI[k] = RefSigsI[k-1] + Nt;
                }

        // Указатель на массив InLLR
            InLLR = mxGetPr(prhs[8]);

    // Выделение памяти под выходную переменную - массив OutLLR. Это
	// выделение памяти надо делать особой функцией, чтобы MATLAB смог
    // получить результат.
        plhs[0] = mxCreateDoubleMatrix(NNSymb*log2M, 1, mxREAL);
        OutLLR = mxGetPr(plhs[0]);
            
    // Выделение памяти для переменных R, Alpha, Beta, Lambda, InSymbsPr
    // Если на выход требуется один аргумент, то все эти переменные
    // внутренние, иначе они объявляются как выходные переменные
    // При выделении памяти через double содержимое выделенной памяти может
    // быть любым, при выделении через mxCreateDoubleMatrix производится
    // заполнение нулями! Алгоритм построен так, что R надо обнулять заново
    // для каждого символа; Alpha, Beta - нужно инициализировать только
    // одно значение и только один раз; Lambda - не нужно инициализировать;
    // InSymbsPr - значения инициализируются внутри LLR2SymbsPr
        if (nlhs == 1) {
            // Указатель на массив R (размер ML x NL1)
                R = new double*[NL1];
                for (k = 0; k < NL1; k++)
                    R[k] = new double[ML];

            // Инициализация указателей на массивы Alpha, Beta и Lambda
            // (размер ML1 x NL), включая инициализацию значений Alpha,
            // Beta
                Alpha = new double*[NL];
                for (k = 0; k < NL; k++)
                    Alpha[k] = new double[ML1];
                *(Alpha[0] + 0) = 1; // По-разному в True Max MaxStar

                Beta = new double*[NL];
                for (k = 0; k < NL; k++)
                    Beta[k] = new double[ML1];
                *(Beta[NL1] + 0) = 1; // По-разному в True Max MaxStar

                Lambda = new double*[NL];
                for (k = 0; k < NL; k++)
                    Lambda[k] = new double[ML1];

            // Указатель на массив InSymbsPr (размер N*NSymb x M)
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
            *(Alpha[0] + 0) = 1; // По-разному в True Max MaxStar

            plhs[3] = mxCreateDoubleMatrix(ML1, NL, mxREAL);
            Beta = new double*[NL];
            Beta[0] = mxGetPr(plhs[3]);
            for (k = 1; k < NL; k++)
                Beta[k] = Beta[k-1] + ML1;
            *(Beta[NL1] + 0) = 1; // По-разному в True Max MaxStar

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

    // Первый этап обработки - преобразование входных llr в вероятности
    // символов для всего кадра
        LLR2SymbsPr(M, log2M, NNSymb, InLLR, InSymbsPr);

    // Цикл по количеству NOFDM-символов
        SpectrsShift   = 0; // сдвиг внутри массивов Spectrs,
        InSymbsPrShift = 0; // InSymbsPr и
        BitNumberShift = 0; // OutLLR
        for (k = 0; k < NSymb; k++) {
            // Расчёт R - условных вероятностей перехода опорных сигналов в
            // принятые сигналы
                SetR2Zero(NL1, ML, R);
                CalculateR(L, N, Nt, ML, ML1, log2M, SpectrsShift, Var,
                        isSigComplex, SpectrsR, SpectrsI, RefSigsR,
                        RefSigsI, R);

            // Расчёт массива Alpha
                CalculateAlpha(L, N, M, ML1, ML2, log2M, InSymbsPrShift,
                        Alpha, InSymbsPr, R);

            // Расчёт массива Beta
                CalculateBeta (L, N, M, ML1, ML2, log2M, InSymbsPrShift,
                        Beta, InSymbsPr, R);

            // Расчёт llr бит
                CalculateBitsLLR(L, N, M, ML, ML1, log2M, BitNumberShift,
                        OutLLR, Alpha, Beta, Lambda);

            // Обновим значение сдвига по NOFDM-символам
                SpectrsShift   += NL1; // сдвиг внутри массивов Spectrs,
                InSymbsPrShift += N;   // InSymbsPr и
                BitNumberShift += N*log2M; // OutLLR
        }

    // Всю выделенную в функции mexFunction память (кроме выходных
    // переменных) надо освободить
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
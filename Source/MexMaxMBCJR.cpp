#include "mex.h"
#include <memory>
#include "../Include/Types.h"
#include "../Include/AlphaBetaMaxMN.h"
#include "../Include/LLR_MN.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Внимание! В функции нет никакой защиты от дурака, поэтому при
    // "вылетании" функции прежде всего надо проверять типы и размеры
    // входных переменных!
    //
    // Пример вызова mex-функции из MATLAB:
    //
    //   OutLLR = MexMaxBCJR(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
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
// Указатель на массив OutLLR
    double *OutLLR;

// Внутренние переменные
// Часто используемые значения: log2M = log2(M), ML = M^L;
// ML1 = M^(L-1); ML2 = M^(L-2); NL = N+L; NL1 = N+L-1,
// NNSymb = N*NSymb
    int log2M, ML, ML1, ML2, NL, NL1, NNSymb;
// Счётчики
    int k, SpectrsShift, InSymbsLogPrShift, BitNumberShift;
// Флаг, указывающий на то комплексный сигнал или нет
    bool isSigComplex;
           
	int MN, MNtype;
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
	// Algorithm
	MN = ( int )( mxGetScalar(prhs[ 9 ]) );
		 
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
	// Указатели на массив Spectrs
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

// Указатель на массив InSymbsPr
	InLLR = mxGetPr(prhs[8]);

	AlgorithmType   aType = MaxMBCJR;

	std::unique_ptr<AlphaBetaMaxMN> alphabeta( new	AlphaBetaMaxMN( N, M, L, aType, MN, Nt, Var ) );
	std::unique_ptr<MnLLR>		 mnllr ( new	MnLLR( N, M, L, NSymb, aType, MN ) );

    // Выделение памяти под выходную переменную - массив OutLLR. Это
	// выделение памяти надо делать особой функцией, чтобы MATLAB смог
    // получить результат.
    plhs[0] = mxCreateDoubleMatrix(NNSymb*log2M, 1, mxREAL);
    OutLLR = mxGetPr(plhs[0]);

    // Выделение памяти для переменных R, Alpha, Beta, Lambda, InSymbsLogPr
    // Если на выход требуется один аргумент, то все эти переменные
    // внутренние, иначе они объявляются как выходные переменные
    // При выделении памяти через double содержимое выделенной памяти может
    // быть любым, при выделении через mxCreateDoubleMatrix производится
    // заполнение нулями! Алгоритм построен так, что R надо обнулять заново
    // для каждого символа; Alpha, Beta - нужно инициализировать только
    // одно значение и только один раз; Lambda - не нужно инициализировать;
    // InSymbsLogPr - значения инициализируются внутри LLR2SymbsLogPr


    // Первый этап обработки - преобразование входных llr в логарифмы
    // вероятностей символов для всего кадра
	 
	mnllr->LLR2SymbsLogPr(InLLR);
    DMatrix&	SymbLogPr = mnllr->GetInSymbLogLLR( );
	 
    SpectrsShift      = 0; // сдвиг внутри массивов Spectrs,
    InSymbsLogPrShift = 0; // InSymbsLogPr и
    BitNumberShift    = 0; // OutLLR
    for (k = 0; k < NSymb; ++k ) {     // Цикл по количеству NOFDM-символов
 
		// Расчёт Alpha
		alphabeta->CalculateAlpha(SymbLogPr, SpectrsR, RefSigsR, InSymbsLogPrShift, SpectrsShift);
		// Расчёт Beta
		alphabeta->CalculateBeta(SymbLogPr, SpectrsR, RefSigsR, InSymbsLogPrShift, SpectrsShift);

		StMatrix& Alpha = alphabeta->GetAlpha();
		StMatrix& Beta = alphabeta->GetBeta();

		// Расчёт llr бит
		mnllr->CalculateBitsMaxLLR( &OutLLR[ BitNumberShift ], Alpha, Beta );

		// Обновим значение сдвига по NOFDM-символам
		SpectrsShift      += NL1; // сдвиг внутри массивов Spectrs
		InSymbsLogPrShift += N;   // InSymbsLogPr и
		BitNumberShift    += N*log2M; // OutLLR*/
    } 
    // Всю выделенную в функции mexFunction память (кроме выходных
    // переменных) надо освободить
 
    delete [] SpectrsR;
    delete [] RefSigsR;
    if (isSigComplex){
        delete [] SpectrsI;
        delete [] RefSigsI;
    }

    return;
}




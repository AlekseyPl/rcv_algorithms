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
// ФУНКЦИЯ-ОБЁРТКА //////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Внимание! В функции нет никакой защиты от дурака, поэтому при
    // "вылетании" функции прежде всего надо проверять типы и размеры
    // входных переменных!
    //
    // Пример вызова mex-функции из MATLAB:
    //
    //   OutLLR = MexViterbi(L, M, N, Nt, NSymb, Var, Spectrs, RefSigs,
    //          InLLR, TBLen);
    //      L     - кол-во компонент, натуральное число;
    //      M     - размер созвездия, натуральное число;
    //      N     - кол-во модуляционных символов в одном NOFDM-символе,
    //              натуральное число;
    //      Nt    - кол-во отсчётов сигнала на один ТИ, натуральное число;
    //      NSymb - кол-во NOFDM-символов в одном кадре, натуральное число;
    //      Var   - дисперсия АБГШ, вещественное число; - НЕ ИСПОЛЬЗУЕТСЯ
    //      Spectrs - (Nt*(N+L-1) x NSymb) массив вещественных /
    //                комплексных спектров - здесь он будет трактоваться
    //                как массив (Nt x (N+L-1)*NSymb);
    //      RefSigs - (Nt х (M^L)*(2*L-1)) массив вещественных /
    //                комплексных опорных сигналов;
    //      InLLR   - (N*NSymb*log2(M) х 1) массив вещественных значений
    //                априорных llr бит; - НЕ ИСПОЛЬЗУЕТСЯ, можно подавать
    //                всё, что угодно;
    //      TBLen   - возможные значения от 1 до (N+L-1).
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
		int NSymb;

	// Указатели на вещественную (и мнимую) часть
		DMatrix SpectrsR;
		DMatrix SpectrsI;
		DMatrix RefSigsR;
		DMatrix RefSigsI;


	// Значение длины памяти демодулятора, т.е. количества тактовых
	// интервалов сигнала на основании которых получается значение
	// первого демодулированного символа
		int TBLen;

// Выходные переменные
		int k, SpectrsShift, BitNumShift;

		SignalParams params;
    // Получение значений входных переменных
	// Количество компонент
		params.L  = static_cast<int>(mxGetScalar(prhs[L_idx]));
	// Размер созвездия
		params.M  = static_cast<int>(mxGetScalar(prhs[M_idx]));
	// Кол-во модуляционных символов в одном NOFDM-символе
		params.N  = static_cast<int>(mxGetScalar(prhs[N_idx]));
	// Кол-во отсчётов сигнала на один ТИ
		params.Nt = static_cast<int>(mxGetScalar(prhs[Nt_idx]));
	// Кол-во NOFDM-символов в одном кадре
		NSymb = static_cast<int>(mxGetScalar(prhs[Nsymb_idx]));
	// Дисперсия шума - НЕ ИСПОЛЬЗУЕТСЯ!
	//    Var = (double)(mxGetScalar(prhs[5]));
	// TBLen
		TBLen = static_cast<int>(mxGetScalar(prhs[TBLen_idx]));

    // Расчёт часто используемых значений - для упрощения кода    
         
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

	// Флаг, указывающий на то комплексный сигнал или нет
		bool isSigComplex = mxIsComplex(prhs[Spectrs_idx]);
		params.type = isSigComplex ? DataType::complex : DataType::real;
		params.constType = ( params.M > 2 ) ? ConstellationType::m_ary : ConstellationType::binary;


// Продолжение получения значений входных переменных

		SpectrsR.resize(params.NL1*NSymb);
		auto ptrR = mxGetPr(prhs[Spectrs_idx]);
		for( int i =0; i < SpectrsR.size(); ++i ) {
			SpectrsR.resize(params.Nt);
			for( auto& v : SpectrsR[ i ] ) v = *ptrR++;
		}


		if (isSigComplex) {	// Мнимая часть
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

		if (isSigComplex) {		// Мнимая часть
			RefSigsI.resize(params.ML2L1);
			auto ptrRefI = mxGetPi(prhs[RefSigs_idx]);
			for( int i =0; i < params.ML2L1; ++i ) {
				RefSigsI.resize(params.Nt);
				for( auto& v : RefSigsI[ i ] ) v = *ptrRefI++;
			}
		}


    // Выделение памяти под выходную переменную - массив OutLLR. Это
	// выделение памяти надо делать особой функцией, чтобы MATLAB смог
    // получить результат.
		plhs[0] = mxCreateDoubleMatrix(NNSymb*params.log2M, 1, mxREAL);
		double *OutLLR = mxGetPr(plhs[0]);
		std::vector<double> vecOut;
		vecOut.reserve(NNSymb*params.log2M);

		shared_ptr< ViterbiDemodulator > vitD = make_shared< ViterbiDemodulator >( params, TBLen );

    // Цикл по количеству NOFDM-символов
		SpectrsShift = 0; // сдвиг внутри массивов Spectrs и
		BitNumShift  = 0; // OutLLR
		for (k = 0; k < NSymb; k++) {
			// Вызов основной C-функции
				vitD->process( RefSigsR, RefSigsI, SpectrsR, SpectrsI, vecOut, SpectrsShift, BitNumShift);
			// Обновим значение сдвига по NOFDM-символам
				SpectrsShift += params.NL1; // сдвиг внутри массивов Spectrs и
				BitNumShift  += params.N*params.log2M; // OutLLR
		}
//		memcpy(OutLLR, &vecOut[0], vecOut.size() * sizeof(double));

//		for( auto&i : vecOut) 			*OutLLR++ = i;
//		for( int i = 0;)
        
    return;
}

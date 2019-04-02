#include "mex.h"  
#include "../Include/Types.h"
#include "../Include/AlphaBetaTrueMN.h"
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
    //          InLLR, MN, MNType );
    //
    //   Входные переменные:
    //      L       - кол-во компонент, натуральное число;
    //      M       - размер созвездия, натуральное число;
    //      N       - кол-во модуляционных символов в одном NOFDM-символе,	 натуральное число;
    //      Nt      - кол-во отсчётов сигнала на один ТИ, натуральное число;
    //      NSymb   - кол-во NOFDM-символов в одном кадре, натуральное число;
    //      Var     - дисперсия АБГШ, вещественное число;
    //      Spectrs - (Nt*(N+L-1) x NSymb) массив вещественных /
    //                комплексных спектров - здесь он будет трактоваться
    //                как массив (Nt x (N+L-1)*NSymb);
    //      RefSigs - (Nt х (M^L)*(2*L-1)) массив вещественных /
    //                комплексных опорных сигналов;
    //      InLLR   - (N*NSymb*log2(M) х 1) массив вещественных значений
    //                априорных llr бит;
    //      MN      - Количество выживших путей, натуральное число;
	// 		MNType  - Тип вычисления метрик обратного прохода: лучшее значение/по метрикам прямого прохода
    //
    //   Выходные переменные:
    //      OutLLR  - (N*NSymb*log2(M) х 1) массив вещественных значений
    //                апостериорных llr бит.
    //
    // Внимание к матричным переменным! При передаче двумерного массива,
    // т.е. матрицы, из Matlab в C++ данные записываются в одномерный
    // массив, причём порядок считывания из матрицы - по столбцам!

 
	// целочисленные переменные-скаляры
	int L, M, N, Nt, NSymb;
	int MN;
	// вещественные переменные-скаляры
	double Var;
    // Указатели на вещественную (и мнимую) часть(и) переменных-массивов
	double **SpectrsR, **SpectrsI, **RefSigsR, **RefSigsI,	*InLLR;
	// Указатель на массив OutLLR

	double *OutLLR;
	// Р’РЅСѓС‚СЂРµРЅРЅРёРµ РїРµСЂРµРјРµРЅРЅС‹Рµ
	// Р§Р°СЃС‚Рѕ РёСЃРїРѕР»СЊР·СѓРµРјС‹Рµ Р·РЅР°С‡РµРЅРёСЏ: log2M = log2(M), ML = M^L;
	// ML1 = M^(L-1); ML2 = M^(L-2); NL = N+L; NL1 = N+L-1,
	// NNSymb = N*NSymb
	int log2M, ML, ML1, ML2, NL, NL1, NNSymb;
	// РЎС‡С‘С‚С‡РёРєРё
	int SpectrsShift, InSymbsPrShift, BitNumberShift;
	 // Р¤Р»Р°Рі, СѓРєР°Р·С‹РІР°СЋС‰РёР№ РЅР° С‚Рѕ РєРѕРјРїР»РµРєСЃРЅС‹Р№ СЃРёРіРЅР°Р» РёР»Рё РЅРµС‚
	bool isSigComplex;


// Получение значений входных переменных
	// Количество компонент
	L  = (int)(mxGetScalar(prhs[0]));
	// Р Р°Р·РјРµСЂ СЃРѕР·РІРµР·РґРёСЏ
	M  = (int)(mxGetScalar(prhs[1]));
	// РљРѕР»-РІРѕ РјРѕРґСѓР»СЏС†РёРѕРЅРЅС‹С… СЃРёРјРІРѕР»РѕРІ РІ РѕРґРЅРѕРј NOFDM-СЃРёРјРІРѕР»Рµ
	N  = (int)(mxGetScalar(prhs[2]));
	// РљРѕР»-РІРѕ РѕС‚СЃС‡С‘С‚РѕРІ СЃРёРіРЅР°Р»Р° РЅР° РѕРґРёРЅ РўР�
	Nt = (int)(mxGetScalar(prhs[3]));
	// РљРѕР»-РІРѕ NOFDM-СЃРёРјРІРѕР»РѕРІ РІ РѕРґРЅРѕРј РєР°РґСЂРµ
	NSymb = (int)(mxGetScalar(prhs[4]));
	// Р”РёСЃРїРµСЂСЃРёСЏ С€СѓРјР°
	Var = (double)(mxGetScalar(prhs[5]));
	// Р¤Р»Р°Рі, СѓРєР°Р·С‹РІР°СЋС‰РёР№ РЅР° С‚Рѕ РєРѕРјРїР»РµРєСЃРЅС‹Р№ СЃРёРіРЅР°Р» РёР»Рё РЅРµС‚
	isSigComplex = mxIsComplex(prhs[6]);
	// Algorithm
	MN = (int)(mxGetScalar(prhs[9]));
 
	// Р Р°СЃС‡С‘С‚ С‡Р°СЃС‚Рѕ РёСЃРїРѕР»СЊР·СѓРµРјС‹С… Р·РЅР°С‡РµРЅРёР№ - РґР»СЏ СѓРїСЂРѕС‰РµРЅРёСЏ РєРѕРґР°
	// log2M = log2(M)
	log2M = 1;
	int k = 2;
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
	// -1/(2*Var) - РјРЅРѕР¶РёС‚РµР»СЊ РїРѕРєР°Р·Р°С‚РµР»СЏ СЌРєСЃРїРѕРЅРµРЅС‚С‹
	Var = -1/(2*Var);
    AlgorithmType   aType = TrueMBCJR;
	    
	AlphaBetaTrueMN* alphabeta = new	AlphaBetaTrueMN( N, M, L, aType, MN, Nt, Var );
	MnLLR*	 	 mnLlr 	  = new		MnLLR( N, M, L, NSymb, aType, MN );
 	// РџСЂРѕРґРѕР»Р¶РµРЅРёРµ РїРѕР»СѓС‡РµРЅРёСЏ Р·РЅР°С‡РµРЅРёР№ РІС…РѕРґРЅС‹С… РїРµСЂРµРјРµРЅРЅС‹С…
	// РЈРєР°Р·Р°С‚РµР»Рё РЅР° РјР°СЃСЃРёРІ Spectr
	// Р’РµС‰РµСЃС‚РІРµРЅРЅР°СЏ С‡Р°СЃС‚СЊ
	SpectrsR = new double*[NL1*NSymb];
	SpectrsR[0] = mxGetPr(prhs[6]);
	for ( int  k = 1; k < NL1*NSymb; k++)
		SpectrsR[k] = SpectrsR[k-1] + Nt;
	// РњРЅРёРјР°СЏ С‡Р°СЃС‚СЊ
	if (isSigComplex) {
		SpectrsI = new double*[NL1*NSymb];
		SpectrsI[0] = mxGetPi(prhs[6]);
		for ( int k = 1; k < NL1*NSymb; k++)
			SpectrsI[k] = SpectrsI[k-1] + Nt;
	}
            
	// РЈРєР°Р·Р°С‚РµР»Рё РЅР° РјР°СЃСЃРёРІ RefSigs
	// Р’РµС‰РµСЃС‚РІРµРЅРЅР°СЏ С‡Р°СЃС‚СЊ
	RefSigsR = new double*[ML*(2*L-1)];
	RefSigsR[0] = mxGetPr(prhs[7]);
	for ( int k = 1; k < ML*(2*L-1); k++)
		RefSigsR[k] = RefSigsR[k-1] + Nt;
	// РњРЅРёРјР°СЏ С‡Р°СЃС‚СЊ
	if (isSigComplex) {                
		RefSigsI = new double*[ML*(2*L-1)];
		RefSigsI[0] = mxGetPi(prhs[7]);
		for ( int  k = 1; k < ML*(2*L-1); k++)
			RefSigsI[k] = RefSigsI[k-1] + Nt;
	}
	 
	// РЈРєР°Р·Р°С‚РµР»СЊ РЅР° РјР°СЃСЃРёРІ InLLR
	InLLR = mxGetPr(prhs[8]);

   // Р’С‹РґРµР»РµРЅРёРµ РїР°РјСЏС‚Рё РїРѕРґ РІС‹С…РѕРґРЅСѓСЋ РїРµСЂРµРјРµРЅРЅСѓСЋ - РјР°СЃСЃРёРІ OutLLR. Р­С‚Рѕ
	// РІС‹РґРµР»РµРЅРёРµ РїР°РјСЏС‚Рё РЅР°РґРѕ РґРµР»Р°С‚СЊ РѕСЃРѕР±РѕР№ С„СѓРЅРєС†РёРµР№, С‡С‚РѕР±С‹ MATLAB СЃРјРѕРі
    // РїРѕР»СѓС‡РёС‚СЊ СЂРµР·СѓР»СЊС‚Р°С‚.
	plhs[0] = mxCreateDoubleMatrix(NNSymb*log2M, 1, mxREAL);
	OutLLR = mxGetPr(plhs[0]);

	 // Первый этап обработки - преобразование входных llr в ероятности символов для всего кадра
	mnLlr->LLR2SymbsPr( InLLR  );
	DMatrix& InSymbsPr = mnLlr->GetInSymbLLR( );
    // Цикл по количеству NOFDM-символов
 
    SpectrsShift   = 0; // сдвиг внутри массивов Spectrs,
    InSymbsPrShift = 0; // InSymbsLogPr и
    BitNumberShift = 0; // OutLLR

	for ( int k = 0; k < NSymb; ++k ) {
 
	// Расчёт Alpha
		alphabeta->CalculateAlpha( InSymbsPr, SpectrsR, RefSigsR, InSymbsPrShift, SpectrsShift);
	// Расчёт Beta
		alphabeta->CalculateBeta ( InSymbsPr, SpectrsR, RefSigsR, InSymbsPrShift, SpectrsShift);

		StMatrix& Alpha = alphabeta->GetAlpha( );
		StMatrix& Beta = alphabeta->GetBeta();
 
   // Расчёт llr бит
		mnLlr->CalculateBitsLLR( &OutLLR[ BitNumberShift ], Alpha, Beta );
	 	StMatrix& Lambda = mnLlr->GetLambda();
	 
   // Обновим значение сдвига по NOFDM-символам
		SpectrsShift   += NL1; // сдвиг внутри массивов Spectrs,
		InSymbsPrShift += N;   // InSymbsLogPr и
		BitNumberShift += N*log2M; // OutLLR
	}
	
    // Всю выделенную в функции mexFunction память (кроме выходных
    // переменных) надо освободить
 

	delete [] SpectrsR;
	delete [] RefSigsR;

	if ( isSigComplex ) {
		delete [] SpectrsI;
		delete [] RefSigsI;
	}

    if( alphabeta ) delete alphabeta;
	if( mnLlr )		delete mnLlr;


    return;
}

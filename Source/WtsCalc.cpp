#include "../Include/WtsCalc.h"
#include <math.h>
#include <iostream>
using namespace std;
  
WtsCalc::WtsCalc( const SignalParams& sp, bool debug ):
	ML(sp.ML), Nt( sp.Nt ), L( sp.L ), NSyms(sp.N), M( sp.M ), type( sp.type ),
	ML1(sp.ML1), ML2( sp.ML2 ), log2M( sp.log2M ), debug( debug )
{
	wts.resize( ML );
}

WtsCalc::~WtsCalc()
{

}

// Функция вычисления весов выживших путей и значений предыдущих
// состояний. Особенность: передаётся указатель не на весь массив
// PrevSts, а только на текущий столбец PrevStsCurCol
void WtsCalc::CalculatePrevStAndCumWts( DVector& NewCumWts, DVector& OldCumWts, IVector& PrevStsCurCol)
{
// Для каждого возможного состояния (CurSt) определим лучший путь
// приходящий в него. Первые два цикла (по n и по m) нужны для
// того, чтобы PrevSt менялся особым образом с изменением CurSt.
// Особенности изменения RefSigNum и PrevSt при изменении CurSt от
// 0 до M-1 можно увидеть из "решётки" сигнала.
	int CurSt = 0;
	int PrevSt;
	int RefSigNum = 0;
	double Buf = 0;
	double Min;
	std::fill( NewCumWts.begin(), NewCumWts.end(), MInf );
	std::fill( PrevStsCurCol.begin(), PrevStsCurCol.end(), -1 );

	for ( int n = 0; n < M; ++n ) {
		PrevSt = 0;
		for ( int m = 0; m < ML2; ++m ) {
			Min = MInf;
			// Цикл по состояниям, из которых можно попасть в текущее
			for (int p = 0; p < M; ++p) {
				// Вычислим метрику пути из PrevSt в CurSt
				Buf = OldCumWts[ PrevSt ] + wts[ RefSigNum ];
				// Если она минимальная, то надо её сохранить
				if (Buf < Min) {
					Min = Buf;
					NewCumWts[ CurSt ] = Buf;
					PrevStsCurCol[ CurSt ] = PrevSt;
				}
				// Подготовка к следующей итерации
				++PrevSt;
				++RefSigNum;
			}
			++CurSt;
		}
	}
}
 

void WtsCalc::CalculateWts(const DMatrix& RefSigsR, const DMatrix& RefSigsI,
						   const DVector& SigR, const DVector& SigI, int k )
{
	if (type == DataType::real)		CalculateReal(RefSigsR, SigR, k);
	else		                    CalculateComplex(RefSigsR, RefSigsI, SigR, SigI, k);
}
// Функция вычисления значений весов всех переходов, т.е. расстояний от
// принятого на текущем тактовом интервале сигнала до опорных сигналов.
// Выбор функции между WtsCalculationCAnyNt, WtsCalculationCNtOne,
// WtsCalculationRAnyNt и WtsCalculationRNtOne определяется
// комплексностью сигнала и значением Nt
void WtsCalc::CalculateComplex(const DMatrix& RefSigsR, const DMatrix& RefSigsI,
							   const DVector& SigR, const DVector& SigI, int k)
{
	// Определим сдвиг в массивах RefSigs
	if (k > L - 2) 			RefSigShift = ( k < NSyms ) ? L - 1 : L + k - NSyms;
	else					RefSigShift = k;

	RefSigShift *= ML;

	// Номер текущего тактового интервала - k1
	// Для каждого опорного сигнала n вычисляется евклидово расстояние 
	// до принятого на текущем, т.е. k1-м, тактовом интервале сигнала
	fill(wts.begin(), wts.end(), 0);
	for (auto wtsIt = wts.begin(); wtsIt < wts.end(); wtsIt++) {

		BufRefSigShift = RefSigShift + distance(wts.begin(), wtsIt);
		for (int m = 0; m < Nt; ++m) {
			auto re = RefSigsR[ BufRefSigShift ][ m ] - SigR[ m ];
			auto im = RefSigsI[ BufRefSigShift ][ m ] - SigI[ m ];
			*wtsIt += re * re + im * im;
		}
	}
}

void WtsCalc::CalculateReal(const DMatrix& RefSigsR, const DVector& SigR, int k )
{
	// Определим сдвиг в массивах RefSigs
	if (k > L - 2) 			RefSigShift = ( k < NSyms ) ? L - 1 : L + k - NSyms;
	else					RefSigShift = k;
	RefSigShift *= ML;
	// Номер текущего тактового интервала - k1
	// Для каждого опорного сигнала n вычисляется евклидово расстояние 
	// до принятого на текущем, т.е. k1-м, тактовом интервале сигнала
	fill(wts.begin(), wts.end(), 0);
	BufRefSigShift = RefSigShift;
	for (auto wtsIt = wts.begin(); wtsIt < wts.end(); ++wtsIt) {
		for (int n = 0; n < Nt; ++n ) {
			auto re = RefSigsR[ BufRefSigShift ][ n ] - SigR[ n ];
			*wtsIt += re * re;
		}
		++BufRefSigShift;
	}
}

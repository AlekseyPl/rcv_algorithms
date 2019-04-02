#include "../Include/Sym2Bit.h"
#include <math.h>

// Функция вычисления и записи значения оценки бит при известном
// значении оценки символа. При M = 2 задача вырожденная и используется
// функция Sym2Bits2, в противном случае используется Sym2BitsM.
// Вместо самого BitNum передаётся указатель на него, чтобы изменения
// BitNum в функции Sym2Bits приводили к изменению BitNum и в основной
// функции.

Sym2Bit::Sym2Bit( const SignalParams& sp ):
	log2M( sp.log2M ), constellation( sp.constType )
{

}

Sym2Bit::~Sym2Bit()
{

}

void Sym2Bit::convert(DVector& OutLLR, const int& OutSym )
{
	if (constellation == ConstellationType::binary) binaryConvert(OutLLR, OutSym);
	else											maryConvert(OutLLR, OutSym);
}

void Sym2Bit::binaryConvert( DVector& OutLLR, const int& OutSym )
{
	// Так было - значение бита 0/1
	// *(OutLLR + *BitNum) = (double)(OutSym);
	OutLLR.push_back( 1000 - 2000 * static_cast<double>( OutSym ) );
}

void Sym2Bit::maryConvert( DVector& OutLLR,const int& OutSym )
{
 	// Для случая, когда старший бит идёт первым
	auto Mult = ( 1 << ( log2M - 1 ) );
	for ( int n = log2M - 1; n >= 0; n--) {
		// Так было - значение бита 0/1
		// *(OutLLR + *BitNum) = (double)((OutSym & Mult) >> n);
		// Так стало - значение llr бита 1000/-1000
		OutLLR.push_back(1000 -	2000 * static_cast<double>( ( OutSym & Mult ) >> n ) );
		Mult >>= 1;
	}
}

#if 0 // А тут старое
// Для случая, когда младший бит идёт первым
//         for (n = 0; n < log2M; n++){
//             // Так было - значение бита 0/1
//             // *(OutLLR + *BitNum) = (double)(OutSym & 1);
//             // Так стало - значение llr бита 1000/-1000
//             *(OutLLR + *BitNum) = 1000 - 2000*((double)(OutSym & 1));
//             OutSym >>= 1;
//             (*BitNum)++;
//         }
#endif

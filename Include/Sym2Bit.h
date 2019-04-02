#/*
 * Sym2Bit.h
 *
 *  Created on: Nov 9, 2018
 *      Author: aplotnikov
 */

#ifndef INCLUDE_SYM2BIT_H_
#define INCLUDE_SYM2BIT_H_

#include "Types.h"
#include <algorithm>
#include <vector>
#include <memory>

class Sym2Bit {
public: 

	Sym2Bit( const SignalParams& sp );
	~Sym2Bit();

	void convert(DVector& OutLLR, const int& OutSym);
private:

	int				log2M;
	ConstellationType	constellation;

	void binaryConvert(DVector& OutLLR, const int& OutSym);
	void maryConvert(DVector& OutLLR, const int& OutSym);
};


#endif

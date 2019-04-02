#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <vector>
#include <float.h>

enum class DataType {
	real = 0,
	complex = 1
};

enum class ConstellationType {
	binary = 0,
	m_ary  = 1

};

struct SignalParams 
{ // ���������, ���������� � ���� ��������� ������� � �������������� �����������
  // ��� ����������� ����������

	int M; // ������ ���������
	int L; // ������������ �������� � �������� ����������
	int Nt;// ���������� �������� �� ��
	int N; // ���������� ������������� ��������



	int log2M;
	int ML;   // M^L
	int ML1;  // M^(L-1)
	int ML2;  // M^(L-2)
	int ML2L1; // M^L * (2*L-1)

	int NL;    // N + L
	int NL1;   // N + L - 1
	int NBits;
	DataType type;					// �����������/������������ ������
	ConstellationType constType;	// ��� ���������(�������� ��������� ����������)

};


enum 	AlgorithmType {
	TrueBCJR  = 0,
	MaxBCJR   = 1,
	TrueMBCJR = 2,
	MaxMBCJR  = 3
};
 
struct ValueState {
    unsigned		state;
    double			value;
    ValueState( unsigned st = 0, double val = 0 ):
    	state( st ), value( val ) { };
};

#ifdef BCJR
const double MInf = -DBL_MAX;
#else
const double MInf = DBL_MAX;
#endif

typedef std::vector< ValueState > 		StVector;
typedef std::vector< StVector > 		StMatrix;

typedef std::vector< double > 			DVector;
typedef	std::vector< DVector > 			DMatrix;

typedef std::vector< int >				IVector;
typedef std::vector< IVector >			IMatrix;


#endif

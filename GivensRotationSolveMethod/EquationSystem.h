#pragma once

#define TEST_UNIT_VECTOR false
#define DEBUG_MODE false

#include <math.h>
#include <complex>
#include <iostream>

using namespace std::complex_literals;


class EquationSystem
{
private:
	int m_SystemSize;
	
	std::complex<double> ** m_A; // Left part of system
	std::complex<double> * m_F;  // Right part of system
	
	std::complex<double> * X_Generated;
	std::complex<double> * X_Found;

	void AllocateCells( int );
	
	void GetRotationMatrix( std::complex<double> &C, std::complex<double> &S, int i, int j );
	
	void MultiplySystem( int i, int j, std::complex<double> C, std::complex<double> S );

	double GetRandomDouble( double dMin, double dMax );
	int GetRandomInt( int iMin, int iMax );

	std::complex<double> * MatrixByVector( std::complex<double> ** Matrix, std::complex<double> * Vector );

public:
	EquationSystem( int, double MinRandom = -10, double MaxRandom = 10 );

	void Solve();

	void PrintSystem( int Precision = 6 );

	void PrintVector( std::complex<double> * Vector, int Precision = 6 );
};

EquationSystem::EquationSystem( int Size, double MinRandom, double MaxRandom )
{
	AllocateCells( Size );

	for ( int i = 0; i < m_SystemSize; i++ )
	{
		/* Generating X vector */
		if ( TEST_UNIT_VECTOR )
		{
			X_Generated [ i ].imag( 1 );
			X_Generated [ i ].real( 1 );
		}
		else
		{
			if( !DEBUG_MODE )
			{
				X_Generated [ i ].imag( GetRandomDouble( MinRandom, MaxRandom ) );
				X_Generated [ i ].real( GetRandomDouble( MinRandom, MaxRandom ) );
			}
			else
			{
				X_Generated [ i ].imag( GetRandomInt( MinRandom, MaxRandom ) );
				X_Generated [ i ].real( GetRandomInt( MinRandom, MaxRandom ) );
			}
		}

		/* Generating A matrix */
		for ( int j = 0; j < m_SystemSize; j++ )
		{
			if( !DEBUG_MODE )
			{
				m_A [ i ] [ j ].imag( GetRandomDouble( MinRandom, MaxRandom ) );
				m_A [ i ] [ j ].real( GetRandomDouble( MinRandom, MaxRandom ) );
			}
			else
			{
				m_A [ i ] [ j ].imag( GetRandomInt( MinRandom, MaxRandom ) );
				m_A [ i ] [ j ].real( GetRandomInt( MinRandom, MaxRandom ) );
			}
		}
	}

	/* Solving F vector */
	m_F = MatrixByVector( m_A, X_Generated );

	PrintSystem();
}

void EquationSystem::AllocateCells( int Size )
{
	m_SystemSize = Size;

	m_F = new std::complex<double> [ m_SystemSize ];
	X_Found = new std::complex<double> [ m_SystemSize ];
	X_Generated = new std::complex<double> [ m_SystemSize ];

	m_A = new std::complex<double> * [ m_SystemSize ];

	for ( int i = 0; i < m_SystemSize; i++ )
	{
		m_A [ i ] = new std::complex<double> [ m_SystemSize ];
	}
}

double EquationSystem::GetRandomDouble( double dMin, double dMax )
{
	double d = ( double ) rand() / RAND_MAX;

	double RandomValue = dMin + d * ( dMax - dMin );

	return RandomValue ? RandomValue : 1;
}

int EquationSystem::GetRandomInt( int iMin, int iMax )
{
	int d = ( int ) rand() % ( iMax - iMin );
	int RandomValue = iMin + d;
	return RandomValue ? RandomValue : 1;
}

std::complex<double> * EquationSystem::MatrixByVector( std::complex<double> ** Matrix, std::complex<double> * Vector )
{
	std::complex<double> * temp;
	temp = new std::complex<double> [ m_SystemSize ];

	for ( int i = 0; i < m_SystemSize; i++ )
	{		
		temp [ i ] = 0;

		for ( int j = 0; j < m_SystemSize; j++ )
		{
			temp [ i ] += Matrix [ i ] [ j ] * Vector [ j ];
		}
	}

	return temp;
}

/* Fills params C and S with required elements to allow A * T = 0, where T is rotation matrix */
void EquationSystem::GetRotationMatrix( std::complex<double> &C, std::complex<double> &S, int i, int j )
{
	std::complex<double> buf = sqrt( m_A [ i ] [ i ] * m_A [ i ] [ i ] + m_A [ j ] [ i ] * m_A [ j ] [ i ] );

	C = m_A [ i ] [ i ] / buf;
	S = m_A [ j ] [ i ] / buf;

}

void EquationSystem::MultiplySystem( int i, int j, std::complex<double> C, std::complex<double> S )
{
	/* Handling with right part of the system */
	std::complex<double> Buf_I = m_F [ i ];
	std::complex<double> Buf_J = m_F [ j ];

	m_F [ i ] = C *Buf_I + S * Buf_J;
	m_F [ j ] = C * Buf_J - S * Buf_I;

	/* Handling with left part of the system */
	for ( int k = 0; k < m_SystemSize; k++ )
	{
		Buf_I = m_A [ i ] [ k ];
		Buf_J = m_A [ j ] [ k ];

		m_A [ i ] [ k ] = C * Buf_I + S * Buf_J;
		m_A [ j ] [ k ] = C * Buf_J - S * Buf_I;
	}
}

void EquationSystem::Solve()
{
	for ( int i = 0; i < m_SystemSize - 1; i++ )
	{
		for ( int j = i + 1; j < m_SystemSize; j++ )
		{
			std::complex<double> C;
			std::complex<double> S;

			GetRotationMatrix( C, S, i, j );

			MultiplySystem( i, j, C, S );
		}
	}

	X_Found [ m_SystemSize - 1 ] = m_F [ m_SystemSize - 1 ] / m_A [ m_SystemSize - 1 ] [ m_SystemSize - 1 ];
	for ( int i = m_SystemSize - 2; i >= 0; i-- )
	{
		std::complex<double> temp( 0, 0 );

		for ( int j = i + 1; j < m_SystemSize; j++ )
		{
			temp += X_Found [ j ] * m_A [ i ] [ j ];
		}

		X_Found [ i ] = ( m_F [ i ] - temp ) / m_A [ i ] [ i ];
	}

	std::cout << "Generated X:\n";
	PrintVector( X_Generated );

	std::cout << "Found X:\n";
	PrintVector( X_Found );

}

void EquationSystem::PrintSystem( int RequiredPrecision  )
{
	for ( int i = 0; i < m_SystemSize; i++ )
	{
		for ( int j = 0; j < m_SystemSize; j++ )
		{
			std::cout.precision( RequiredPrecision );
			std::cout.width( RequiredPrecision * 2 + 5 );

			std::cout << m_A [ i ] [ j ];
		}

		std::cout << "\t\t" << m_F [ i ] << "\n\n";
	}

	std::cout << "\n\n\n";
}

void EquationSystem::PrintVector( std::complex<double> * Vector, int Precision /*= 6 */ )
{
	for ( int i = 0; i < m_SystemSize; i++ )
	{
		std::cout.precision( Precision );
		std::cout.width( Precision * 2 + 5 );

		std::cout << Vector [ i ];
	}

	std::cout << "\n\n";
}

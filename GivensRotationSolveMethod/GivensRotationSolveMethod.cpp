#include "stdafx.h"
#include "EquationSystem.h"

int main()
{
	EquationSystem System( 100 );
	
	System.Solve();

	//System.PrintSystem( 4 );

	system( "pause" );

	return 0;
}


#include "stdafx.h"
#include "EquationSystem.h"

int main()
{
	EquationSystem System( 5 );
	
	System.Solve();

	//System.PrintSystem( 4 );

	system( "pause" );

	return 0;
}


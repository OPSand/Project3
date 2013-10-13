// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"

int _tmain(int argc, _TCHAR* argv[])
{
	// flags
	const int DIM = 2;
	const bool FIXED_SUN = true;

	// initialize solar system
	SolarSystem system = SolarSystem(DIM);
	
	CelestialBody sun = CelestialBody("Sun", 100, &system);
	sun.fixed = FIXED_SUN;
	
	CelestialBody earth = CelestialBody("Earth", 1, &system);
	earth.position(0) = -1.0; // x = -1
	earth.position(1) = 5.0; // y = 5
	earth.velocity.fill(1);

	// debug
	system.setForces();
	for( int i = 0; i < system.n(); i++ )
	{
		cout << system.body(i)->name << endl << endl << "Force:" << endl << system.body(i)->force << endl << "Acc:" << endl << system.body(i)->acc() << endl;
	}
	// end debug

	// TODO: iterate and plot coordinates etc.

	getchar(); // pause

	return 0; // exit
}
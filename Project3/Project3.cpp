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
	earth.position.fill(1);
	earth.velocity.fill(1);

	for( int i = 0; i < system.n(); i++ )
	{
		cout << system.body(i)->name() << endl;
	}

	// TODO: iterate and plot coordinates etc.

	getchar(); // pause

	return 0; // exit
}
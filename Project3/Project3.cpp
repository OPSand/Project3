// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include"armadillo"
#include<list>

using namespace arma;
using namespace std; 

int _tmain(int argc, _TCHAR* argv[])
{
	// flags
	const int DIM = 2;
	const bool FIXED_SUN = true;

	// initialize solar system
	SolarSystem system = SolarSystem();
	
	CelestialBody sun = CelestialBody("Sun", DIM, 100);
	sun.fixed = FIXED_SUN;
	system.add(sun);
	
	CelestialBody earth = CelestialBody("Earth", DIM, 1);
	earth.position.fill(1);
	earth.velocity.fill(1);
	system.add(earth);

	// TODO: iterate and plot coordinates etc.

	getchar(); // pause

	return 0;
}
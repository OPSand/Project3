// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"

vec toCartesian2D(double r, double theta)
{
	vec p = vec(2);
	p[0] = r*cos(theta); // x
	p[1] = r*sin(theta); // y
	return p;
}

int _tmain(int argc, _TCHAR* argv[])
{
	// flags
	const int DIM = 2;
	const bool FIXED_SUN = true;

	// masses
	const double M_SUN = 2e30;
	const double M_EARTH = 6e24;
	const double M_JUPITER = 1.9e27;
	const double M_MARS = 6.6e23;
	const double M_VENUS = 4.9e24;
	const double M_SATURN = 5.5e26;
	const double M_MERCURY = 2.4e23;
	const double M_URANUS = 8.8e25;
	const double M_NEPTUNE = 1.03e26;
	const double M_PLUTO = 1.31e22;

	// distances from Sun (AU)
	const double D_SUN = 0.0;
	const double D_EARTH = 1.0;
	const double D_JUPITER = 5.2;
	const double D_MARS = 1.52;
	const double D_VENUS = 0.72;
	const double D_SATURN = 9.54;
	const double D_MERCURY = 0.39;
	const double D_URANUS = 19.19;
	const double D_NEPTUNE = 30.06;
	const double D_PLUTO = 39.53;

	// initialize solar system
	SolarSystem system = SolarSystem(DIM);
	
	CelestialBody sun = CelestialBody("Sun", M_SUN, &system);
	sun.fixed = FIXED_SUN;
	
	CelestialBody earth = CelestialBody("Earth", M_EARTH, &system);
	const double THETA_EARTH = (4.0/6.0) * math::pi();
	earth.position = toCartesian2D(D_EARTH, THETA_EARTH);
	earth.velocity = toCartesian2D(100, (THETA_EARTH + 0.5 * math::pi()));

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
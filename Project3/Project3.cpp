// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"

// converts polar coordinates to cartesian coordinates
vec toCartesian2D(double r, double theta)
{
	vec p = vec(2);
	p[0] = r*cos(theta); // x
	p[1] = r*sin(theta); // y
	return p;
}

// returns an angle in radians that is orthogonal to theta
double orthogonal2D(double theta, bool clockwise = false)
{
	if( clockwise )
	{
		return (theta - 0.5 * math::pi());
	}
	else // counterclockwise
	{
		return (theta + 0.5 * math::pi());
	}
}

// sets a position vector of length d at angle theta, and
// an orthogonal velocity vector of length v
// to help with initialization of celestial bodies
void initial2D(CelestialBody* cb, double d, double v, double theta)
{
	cb->position = toCartesian2D(d, theta);
	cb->velocity = toCartesian2D(v, orthogonal2D(theta)); // counterclockwise
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

	// polar angles (position) in radians
	const double THETA_EARTH = (4.0/6.0) * math::pi();
	const double THETA_JUPITER = 0.0 * math::pi();
	const double THETA_MARS = 0.0 * math::pi();
	const double THETA_VENUS = 0.0 * math::pi();
	const double THETA_SATURN = 0.0 * math::pi();
	const double THETA_MERCURY = 0.0 * math::pi();
	const double THETA_URANUS = 0.0 * math::pi();
	const double THETA_NEPTUNE = 0.0 * math::pi();
	const double THETA_PLUTO = 0.0 * math::pi();

	// initial velocities (absolute value) in AU/s
	const double V_SUN = 0.0;
	const double V_EARTH = 1.0;
	const double V_JUPITER = 1.0;
	const double V_MARS = 1.0;
	const double V_VENUS = 1.0;
	const double V_SATURN = 1.0;
	const double V_MERCURY = 1.0;
	const double V_URANUS = 1.0;
	const double V_NEPTUNE = 1.0;
	const double V_PLUTO = 1.0;

	// initialize solar system
	SolarSystem system = SolarSystem(DIM);
	
	CelestialBody sun = CelestialBody("Sun", M_SUN, &system);
	sun.fixed = FIXED_SUN;
	
	CelestialBody earth = CelestialBody("Earth", M_EARTH, &system);
	initial2D(&earth, D_EARTH, V_EARTH, THETA_EARTH);

	CelestialBody jupiter = CelestialBody("Jupiter", M_JUPITER, &system);
	initial2D(&jupiter, D_JUPITER, V_JUPITER, THETA_JUPITER);

	// debug
	system.setForces();
	for( int i = 0; i < system.n(); i++ )
	{
		cout << system.body(i)->name << endl << endl << "Force:" << endl << system.body(i)->force << endl << "Acc:" << endl << system.body(i)->acc() << endl << "Pos:" << endl << system.body(i)->position << endl << "Vel:" << endl << system.body(i)->velocity << endl;
	}
	// end debug

	// TODO: iterate and plot coordinates etc.

	getchar(); // pause

	return 0; // exit
}
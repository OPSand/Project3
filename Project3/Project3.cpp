// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"
#include <random>

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
		return (theta - 0.5 * cPI);
	}
	else // counterclockwise
	{
		return (theta + 0.5 * cPI);
	}
}

// return a random floating point number in the interval [minIncl, maxExcl).
double randDbl(double minIncl, double maxExcl, minstd_rand* eng)
{
	uniform_real<double> rd(0.0, (2*cPI));
	return rd(*eng);
}

// sets a position vector of length d at angle theta, and
// an orthogonal velocity vector of length v
// to help with initialization of celestial bodies
void initial2D(CelestialBody* cb, double d, double v, minstd_rand* eng, double theta = -1.0)
{
	if( theta == -1.0 ) // no angle given
	{
		// create random angle
		theta = randDbl(0.0, (2*cPI), eng);
	}
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

	// distances from Sun (m)
	const double D_SUN = 0.0 * cAU;
	const double D_EARTH = 1.0 * cAU;
	const double D_JUPITER = 5.2 * cAU;
	const double D_MARS = 1.52 * cAU;
	const double D_VENUS = 0.72 * cAU;
	const double D_SATURN = 9.54 * cAU;
	const double D_MERCURY = 0.39 * cAU;
	const double D_URANUS = 19.19 * cAU;
	const double D_NEPTUNE = 30.06 * cAU;
	const double D_PLUTO = 39.53 * cAU;

	/* polar angles (position) in radians (not used)
	const double THETA_EARTH = 0.0 * cPI;
	const double THETA_JUPITER = 0.0 * cPI;
	const double THETA_MARS = 0.0 * cPI;
	const double THETA_VENUS = 0.0 * cPI;
	const double THETA_SATURN = 0.0 * cPI;
	const double THETA_MERCURY = 0.0 * cPI;
	const double THETA_URANUS = 0.0 * cPI;
	const double THETA_NEPTUNE = 0.0 * cPI;
	const double THETA_PLUTO = 0.0 * cPI; */

	// initial velocities (absolute value) in m/s
	const double V_SUN = 0.0;
	const double V_EARTH = 29.8e3;
	const double V_JUPITER = 13.1e3;
	const double V_MARS = 24.1e3;
	const double V_VENUS = 35.0e3;
	const double V_SATURN = 9.7e3;
	const double V_MERCURY = 47.9e3; // Don't stop me now...
	const double V_URANUS = 6.8e3;
	const double V_NEPTUNE = 5.4e3;
	const double V_PLUTO = 4.7e3;

	// intitialze random number engine
	minstd_rand eng;
	int randomSeed = (int) clock();
	eng.seed(randomSeed);

	// initialize solar system
	SolarSystem system = SolarSystem(DIM);
	
	CelestialBody sun = CelestialBody("Sun", M_SUN, &system, FIXED_SUN);
	
	CelestialBody mercury = CelestialBody("Mercury", M_MERCURY, &system);
	initial2D(&mercury, D_MERCURY, V_MERCURY, &eng);

	CelestialBody venus = CelestialBody("Venus", M_VENUS, &system);
	initial2D(&venus, D_VENUS, V_VENUS, &eng);

	CelestialBody earth = CelestialBody("Earth", M_EARTH, &system);
	initial2D(&earth, D_EARTH, V_EARTH, &eng);

	CelestialBody mars = CelestialBody("Mars", M_MARS, &system);
	initial2D(&mars, D_MARS, V_MARS, &eng);

	CelestialBody jupiter = CelestialBody("Jupiter", M_JUPITER, &system);
	initial2D(&jupiter, D_JUPITER, V_JUPITER, &eng);

	CelestialBody saturn = CelestialBody("Saturn", M_SATURN, &system);
	initial2D(&saturn, D_SATURN, V_SATURN, &eng);

	CelestialBody uranus = CelestialBody("Uranus", M_URANUS, &system);
	initial2D(&uranus, D_URANUS, V_URANUS, &eng);

	CelestialBody neptune = CelestialBody("Neptune", M_NEPTUNE, &system);
	initial2D(&neptune, D_NEPTUNE, V_NEPTUNE, &eng);

	CelestialBody pluto = CelestialBody("Pluto", M_PLUTO, &system);
	initial2D(&pluto, D_PLUTO, V_PLUTO, &eng);

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
// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"
#include <random>

#pragma region Initialization helper methods
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
	*(cb->position) = toCartesian2D(d, theta);
	*(cb->velocity) = toCartesian2D(v, orthogonal2D(theta)); // counterclockwise
}
#pragma endregion

// behold the almighty main method!
int _tmain(int argc, _TCHAR* argv[])
{
#pragma region Flags and Settings
	// dimensions
	const int DIM = 2;

	// time steps
	const int N_STEPS = 300 * 366; // number of steps
	const double STEP = 24 * 60 * 60; // time step length (s)
	const int PLOT_EVERY = 1; // plot every ...th step
	const int N_PLOT = (N_STEPS / PLOT_EVERY); // how many steps we actually plot

	// compiler flags
	const bool ADD_JUPITER = true;
	const bool ADD_ALL = true; // include the other 7 planets
	const bool DEBUG = false; // use for debugging only
	const bool USE_EULER = true; // use Euler-Cromer method (for comparison)
	const bool USE_RK4 = true; // use Runge-Kutta method

	// if set to true, the Sun will never move
	const bool FIXED_SUN = false;

	// mass multiplier for Jupiter (1.0 = normal)
	const double MEGA_JUPITER = 1.0;
#pragma endregion

#pragma region Constants
	// masses
	const double M_SUN = 2e30;
	const double M_EARTH = 6e24;
	const double M_JUPITER = 1.9e27 * MEGA_JUPITER;
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

	/* polar angles (position) in radians (to override random values if needed)
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
#pragma endregion

#pragma region Initialize Solar System
	// intitialze random number engine
	minstd_rand eng;
	int randomSeed = (int) clock();
	eng.seed(randomSeed);

	// initialize solar system
	SolarSystem system = SolarSystem(DIM, N_STEPS, N_PLOT);
	
	CelestialBody* sun = new CelestialBody("Sun", M_SUN, &system, FIXED_SUN);
	
	CelestialBody* earth = new CelestialBody("Earth", M_EARTH, &system);
	initial2D(earth, D_EARTH, V_EARTH, &eng);

	if( ADD_JUPITER ) 
	{
		CelestialBody* jupiter = new CelestialBody("Jupiter", M_JUPITER, &system);
		initial2D(jupiter, D_JUPITER, V_JUPITER, &eng);
	}

	if( ADD_ALL )
	{
		CelestialBody* mercury = new CelestialBody("Mercury", M_MERCURY, &system);
		initial2D(mercury, D_MERCURY, V_MERCURY, &eng);

		CelestialBody* venus = new CelestialBody("Venus", M_VENUS, &system);
		initial2D(venus, D_VENUS, V_VENUS, &eng);
		
		CelestialBody* mars = new CelestialBody("Mars", M_MARS, &system);
		initial2D(mars, D_MARS, V_MARS, &eng);

		CelestialBody* saturn = new CelestialBody("Saturn", M_SATURN, &system);
		initial2D(saturn, D_SATURN, V_SATURN, &eng);

		CelestialBody* uranus = new CelestialBody("Uranus", M_URANUS, &system);
		initial2D(uranus, D_URANUS, V_URANUS, &eng);

		CelestialBody* neptune = new CelestialBody("Neptune", M_NEPTUNE, &system);
		initial2D(neptune, D_NEPTUNE, V_NEPTUNE, &eng);

		CelestialBody* pluto = new CelestialBody("Pluto", M_PLUTO, &system);
		initial2D(pluto, D_PLUTO, V_PLUTO, &eng);
	}

	if( ! FIXED_SUN )
	{
		*(sun->velocity) = (-system.totalMomentum() / sun->mass); // v = p/m;
		//assert( norm(system.totalMomentum(), DIM) == 0.0 ); // check that total momentum is actually 0
	}
#pragma endregion

#pragma region Iterate and create plots

	if( USE_EULER ) // Euler-Cromer method (for comparison)
	{
		SolarSystem euler = system; // make deep copy of solar system
		int n = euler.n(); // number of celestial bodies

		// iterate and plot coordinates
		for( int i = 0; i < N_STEPS; i++ ) // for each time step
		{
			euler.plotCurrentPositions( i % PLOT_EVERY == 0 ); // if we want to plot this step, do it

			// calculate forces/accelerations based on current postions
			euler.setForces();
			
			/* // so, so slow...
			SolarSystem k1 = system; // copy
			k1.diff();
			system += (k1 * STEP); */

			for( int j = 0; j < n; j++ ) // for each celestial body
			{
				CelestialBody* cb = euler.body(j);

				if( ! cb->fixed ) // a fixed celestial body will never move
				{
					// acc -> velocity (Euler-Cromer, for testing only)
					*(cb->velocity) += STEP * cb->acc();

					// velocity -> position (Euler-Cromer, for testing only)
					*(cb->position) += STEP * *(cb->velocity);
				}
			}
		}

		euler.plotDim(0, "Xeuler.dat");
		euler.plotDim(1, "Yeuler.dat");

		cout << "Finished plotting " << N_PLOT << " of " << N_STEPS << " steps (Euler-Cromer)!" << endl;
	}
		
	if( USE_RK4 ) // Runge-Kutta algorithm
	{
		int n = system.n();

		// iterate and plot coordinates
		for( int i = 0; i < N_STEPS; i++ ) // for each time step
		{
			system.plotCurrentPositions( i % PLOT_EVERY == 0 ); // if we want to plot this step, do it

			// matrices for Runge-Kutta values:
			mat k1_v(DIM, n); // k1, velocity
			mat k1_p(DIM, n); // k1, position
			mat k2_v(DIM, n); // k2, velocity (etc.)
			mat k2_p(DIM, n);
			mat k3_v(DIM, n);
			mat k3_p(DIM, n);
			mat k4_v(DIM, n);
			mat k4_p(DIM, n);

			// initial position and velocity for each celestial body
			mat orig_v(DIM, n);
			mat orig_p(DIM, n);

#pragma region k1

			// calculate forces/accelerations based on current postions
			system.setForces();

			for( int j = 0; j < n; j++ ) // for each celestial body
			{
				CelestialBody* cb = system.body(j);

				// store initial position and velocity
				orig_v.col(j) = *(cb->velocity);
				orig_p.col(j) = *(cb->position);

				// save values
				k1_v.col(j) = STEP * cb->acc();
				k1_p.col(j) = STEP * *(cb->velocity);

				// advance to mid-point after k1
				*(cb->velocity) = orig_v.col(j) + 0.5 * k1_v.col(j);
				*(cb->position) = orig_p.col(j) + 0.5 * k1_p.col(j);
			}

#pragma endregion

#pragma region k2

			// calculate forces/accelerations based on current postions
			system.setForces();

			for( int j = 0; j < n; j++ ) // for each celestial body
			{
				CelestialBody* cb = system.body(j);

				// save values
				k2_v.col(j) = STEP * cb->acc();
				k2_p.col(j) = STEP * *(cb->velocity);

				// switch to new mid-point using k2 instead
				*(cb->velocity) = orig_v.col(j) + 0.5 * k2_v.col(j);
				*(cb->position) = orig_p.col(j) + 0.5 * k2_p.col(j);
			}

#pragma endregion

#pragma region k3

			// calculate forces/accelerations based on current postions
			system.setForces();

			for( int j = 0; j < n; j++ ) // for each celestial body
			{
				CelestialBody* cb = system.body(j);

				// save values
				k3_v.col(j) = STEP * cb->acc();
				k3_p.col(j) = STEP * *(cb->velocity);

				// switch to end-point
				*(cb->velocity) = orig_v.col(j) + k3_v.col(j);
				*(cb->position) = orig_p.col(j) + k3_p.col(j);
			}

#pragma endregion

#pragma region k4

			// calculate forces/accelerations based on current postions
			system.setForces();

			for( int j = 0; j < n; j++ ) // for each celestial body
			{
				CelestialBody* cb = system.body(j);

				// save values
				k4_v.col(j) = STEP * cb->acc();
				k4_p.col(j) = STEP * *(cb->velocity);

				// finally, update position and velocity
				*(cb->velocity) = orig_v.col(j) + (1.0/6.0)*(k1_v.col(j) + 2.0 * k2_v.col(j) + 2.0 * k3_v.col(j) + k4_v.col(j));
				*(cb->position) = orig_p.col(j) + (1.0/6.0)*(k1_p.col(j) + 2.0 * k2_p.col(j) + 2.0 * k3_p.col(j) + k4_p.col(j));
			}

#pragma endregion

		}

		system.plotDim(0, "X.dat");
		system.plotDim(1, "Y.dat");

		cout << "Finished plotting " << N_PLOT << " of " << N_STEPS << " steps (Runge-Kutta)!" << endl;
	}
#pragma endregion

	getchar(); // pause

	return 0; // exit
}


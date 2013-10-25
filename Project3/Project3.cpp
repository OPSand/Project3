// Project3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SolarSystem.h"
#include "CelestialBody.h"
#include <random>
#include <assert.h>

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


vec derivative2D (double t, vec& posVel, CelestialBody celestialBody, SolarSystem system)
{
	// Two solutions to match the existing code. 
	vec derivatives = vec(4);
#pragma region Recalculating the forces
	int n = system.n();
	for (int i= 0; i< n; i++)
	{
		if (system.body(i)->name == celestialBody.name);
		else
		{
			double dist = celestialBody.dist(system.body(i));
			vec r = celestialBody.position_diff(system.body(i)); // gives the force the proper direction ... Pourquoi un vecteur ?!!
			derivatives(2) += (cG * celestialBody.mass * system.body(i)->mass / pow(dist, 3.0)) * r(0)*posVel(0); // Newton's law of gravity on x
			derivatives(3) += (cG * celestialBody.mass * system.body(i)->mass / pow(dist, 3.0)) * r(0)*posVel(1); // Newton's law of gravity on y
		}
	}
#pragma endregion

#pragma region Resetting the Forces
	/*
	derivatives[1] = celestialBody.force(t)/celestialBody.mass;*/
#pragma endregion
	
	derivatives[0] = posVel(2); // Derivative of the position in x
	derivatives[1] = posVel(3); // Derivative of the position in y
	return derivatives;
}


// So this is our Runge Kutta 4 algorithm in 2D for now. 
void rk4_2D(int dim, int h, int time,SolarSystem system, CelestialBody currentCelestialBody)
{
#pragma region About the first part
	//So to begin with, we'll only study a system with two points: the sun and the earth.
	// The method to do so is the same as the method to compute this with a lot more elements
	// The only thing changing is the flag set in the beginning of the project3.cpp main
#pragma endregion
	int nbEq = dim * 2; // We have to compute the velocity + position for all our dim.
	vec positions_np1 = vec(2); // on suppose pour l'instant qu'on process en 2D. A arranger après.
	vec f_np1 = vec(nbEq);
	vec k1 = vec(nbEq);
	vec k2 = vec(nbEq);
	vec k3 = vec(nbEq);
	vec k4 = vec(nbEq);
	double halfH = 0.5*h;// We are computing the RG4 with half steps.
	vec f_n = vec(nbEq);
	vec g_n = currentCelestialBody.velocity;
	for (int i= 0 ; i< dim; i++)
	{
		f_n(i) = currentCelestialBody.position(i);
		f_n(dim - 1 + i) = currentCelestialBody.velocity(i);
	}
#pragma region ki Computations
	k1 = derivative2D(time,f_n,currentCelestialBody,system); // k1 
	for (int i=0; i < dim; i++)
	{
		f_n(i) = currentCelestialBody.position(i) + halfH*k1(i);
		f_n(dim - 1 + i) = currentCelestialBody.velocity(i) + halfH*k1(dim - 1 + i);
	}
	k2 = derivative2D(time + halfH,f_n,currentCelestialBody,system ); // k2
	for (int i=0; i < dim; i++)
	{
		f_n(i) = currentCelestialBody.position(i) + halfH*k2(i);
		f_n(dim - 1 + i) = currentCelestialBody.velocity(i) + halfH*k2(dim - 1 + i);
	}
	k3 = derivative2D(time + halfH,f_n,currentCelestialBody,system); // k3
	for (int i=0; i < dim; i++)
	{
		f_n(i) = currentCelestialBody.position(i) + halfH*k3(i);
		f_n(dim - 1 + i) = currentCelestialBody.velocity(i) + halfH*k3(dim - 1 + i);
	}
	k4 = derivative2D(time + halfH,f_n,currentCelestialBody,system); // k4
#pragma endregion

	// And finally, we update the position vector/
	for (int i=0; i<dim;i++)
	{
		f_np1(i) = currentCelestialBody.position(i) + (double)(1/6)*h*(k1(i) + k2(i) + k3(i) + k4(i));
		f_np1(dim - 1 + i) = currentCelestialBody.velocity(i) + double(1/6)*h*(k1(dim - 1 + i) + k2(dim - 1 + i) + k3(dim - 1 + i) + k4(dim - 1 + i));
	}
	// And lastly, we have to update our r in our global struct. Working with a pointer just to be sure that our new params are updated
	system.setForces();
	
	return;
}



int _tmain(int argc, _TCHAR* argv[])
{
	// dimensions
	const int DIM = 2;

	// time steps
	const int N_STEPS = 1000; // number of steps
	const double DELTA_T = 1000; // time step length (s)
	const int PLOT_EVERY = 1; // plot every ...th step

	// flags
	#define ADD_JUPITER	
	#define ADD_ALL	
	//#define DEBUG

	const bool FIXED_SUN = false;

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

	// intitialze random number engine
	minstd_rand eng;
	int randomSeed = (int) clock();
	eng.seed(randomSeed);

	// initialize solar system
	SolarSystem system = SolarSystem(DIM, N_STEPS);
	
	CelestialBody* sun = new CelestialBody("Sun", M_SUN, &system, FIXED_SUN);
	
	CelestialBody* earth = new CelestialBody("Earth", M_EARTH, &system);
	initial2D(earth, D_EARTH, V_EARTH, &eng);

#ifdef ADD_JUPITER
		CelestialBody* jupiter = new CelestialBody("Jupiter", M_JUPITER, &system);
		initial2D(jupiter, D_JUPITER, V_JUPITER, &eng);
#endif

#ifdef ADD_ALL
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

#endif

	if( ! FIXED_SUN )
	{
		sun->velocity = (-system.totalMomentum() / sun->mass); // v = p/m;
		//assert( norm(system.totalMomentum(), DIM) == 0.0 ); // check that total momentum is 0
	}

#ifdef DEBUG
	system.setForces();
	for( int i = 0; i < system.n(); i++ )
	{
		cout << system.body(i)->name << endl << endl << "Force:" << endl << system.body(i)->force << endl << "Acc:" << endl << system.body(i)->acc() << endl << "Pos:" << endl << system.body(i)->position << endl << "Vel:" << endl << system.body(i)->velocity << endl;
	}

	for( int i = 0; i < N_STEPS; i++ )
	{
		vec v(DIM);
		v(0) = i;
		v(1) = 2*i;
		jupiter->position = v;
		if( ! jupiter->plotCurrentPosition() )
		{
			cout << "plot matrix full! :S" << endl;
		}
	}
	cout << jupiter->plot;

	for( int i = 0; i < system.n(); i++ )
	{
		cout << i << ": " << system.body(i)->name << endl; // output planet names
	}

#else // not DEBUG

	// iterate and plot coordinates
	for( int i = 0; i < N_STEPS; i++ ) // for each time step
	{
		// calculate forces/accelerations based on postions
		system.setForces();

		for( int j = 0; j < system.n(); j++ ) // for each celestial body
		{
			CelestialBody* cb = system.body(j);

			if( ! cb->fixed ) // a fixed celestial body will never move
			{
				// acc -> velocity (RK4)
				vec k1 = cb->velocity + DELTA_T * cb->acc();
				// k2 etc.

				// velocity -> position (RK4)

			}

			if( i % PLOT_EVERY == 0 ) // we want to plot this step
			{
				cb->plotCurrentPosition();
			}
		}
	}

	// plot to file for each CB
	for( int j = 0; j < system.n(); j++ )
	{
		system.body(j)->positionToFile(); // saved as "<name>.txt"
	}

	cout << "Finished plotting " << (N_STEPS / PLOT_EVERY) << " of " << N_STEPS << " steps!";

#endif

	getchar(); // pause

	return 0; // exit
}


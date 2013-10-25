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

int _tmain(int argc, _TCHAR* argv[])
{
	// dimensions
	const int DIM = 2;

	// number of time steps
	const int N_STEPS = 1000;

	// flags
	#define ADD_JUPITER	
	#define ADD_ALL	
	#define DEBUG

	const bool FIXED_SUN=false;

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
	
	CelestialBody sun = CelestialBody("Sun", M_SUN, &system, FIXED_SUN);
	
	CelestialBody earth = CelestialBody("Earth", M_EARTH, &system);
	initial2D(&earth, D_EARTH, V_EARTH, &eng);

#ifdef ADD_JUPITER
		CelestialBody jupiter = CelestialBody("Jupiter", M_JUPITER, &system);
		initial2D(&jupiter, D_JUPITER, V_JUPITER, &eng);
#endif

#ifdef ADD_ALL
		CelestialBody mercury = CelestialBody("Mercury", M_MERCURY, &system);
		initial2D(&mercury, D_MERCURY, V_MERCURY, &eng);

		CelestialBody venus = CelestialBody("Venus", M_VENUS, &system);
		initial2D(&venus, D_VENUS, V_VENUS, &eng);
		
		CelestialBody mars = CelestialBody("Mars", M_MARS, &system);
		initial2D(&mars, D_MARS, V_MARS, &eng);

		CelestialBody saturn = CelestialBody("Saturn", M_SATURN, &system);
		initial2D(&saturn, D_SATURN, V_SATURN, &eng);

		//if( true ) // reproduce bug by change of scope
		//{
			CelestialBody uranus = CelestialBody("Uranus", M_URANUS, &system);
			initial2D(&uranus, D_URANUS, V_URANUS, &eng);

			CelestialBody neptune = CelestialBody("Neptune", M_NEPTUNE, &system);
			initial2D(&neptune, D_NEPTUNE, V_NEPTUNE, &eng);

			CelestialBody pluto = CelestialBody("Pluto", M_PLUTO, &system);
			initial2D(&pluto, D_PLUTO, V_PLUTO, &eng);
		//}
#endif

	if( ! FIXED_SUN )
	{
		sun.velocity = (-system.totalMomentum() / sun.mass); // v = p/m;
		//assert( norm(system.totalMomentum(), DIM) == 0.0 ); // check that total momentum is 0
	}

#ifdef DEBUG
	for( int i = 0; i < system.n(); i++ )
	{
		cout << i << ": " << system.body(i)->name << endl; // output planet names
	}

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
		jupiter.position = v;
		if( ! jupiter.plotCurrentPosition() )
		{
			cout << "plot matrix full! :S" << endl;
		}
	}
	cout << jupiter.plot;

#endif

	// TODO: iterate and plot coordinates etc.

	getchar(); // pause

	return 0; // exit
}


// We just compute the derivatives here, and gives it back to the "rkSolver"
vec derivative (double t, double h, vec& position, vec& velocity, CelestialBody celestialBody, SolarSystem* system)
{
	// Two solutions to match the existing code. 
	vec derivatives = vec(2);
#pragma region Recalculating the forces
	int n = system->n();
	for (int i= 0; i< n; i++)
	{
		if (system->body(i)->name == celestialBody.name);
		else
		{
			double dist = celestialBody.dist(system->body(i));
			vec r = celestialBody.position_diff(system->body(i)); // gives the force the proper direction ... Pourquoi un vecteur ?!!
			derivatives(1) += (cG * celestialBody.mass * system->body(i)->mass / pow(dist, 3.0)) * r(0)*position(0); // Newton's law of gravity
		}
	}
#pragma endregion

#pragma region Resetting the Forces
	/*
	derivatives[1] = celestialBody.force(t)/celestialBody.mass;*/
#pragma endregion
	
	derivatives[0] = velocity(t); // Derivative of the position
	return derivatives;
}


// So this is our Runge Kutta 4 algorithm in 2D for now. 
void rk4_2D(int nbSteps, double h, int time, double stepTime, SolarSystem* system, CelestialBody currentCelestialBody)
{
#pragma region About the first part
	//So to begin with, we'll only study a system with two points: the sun and the earth.
	// The method to do so is the same as the method to compute this with a lot more elements
	// The only thing changing is the flag set in the beginning of the project3.cpp main
#pragma endregion
	vec positions_np1 = vec(2); // on suppose pour l'instant qu'on process en 2D. A arranger après.
	vec f_np1 = vec(2);
	vec k1 = vec(2);
	vec k2 = vec(2);
	vec k3 = vec(2);
	vec k4 = vec(2);
	double halfH = 0.5*h;// We are computing the RG4 with half steps.
	vec f_n = currentCelestialBody.position;
	vec g_n = currentCelestialBody.velocity;
	double halfStepTime = stepTime/2;

#pragma region k Computations
	k1 = derivative(time,nbSteps,currentCelestialBody.position,currentCelestialBody.velocity,currentCelestialBody,system); // k1 computation
	for (int i=0; i < nbSteps; i++)
		f_n(i) = currentCelestialBody.position(i) + halfH*k1(i);
	k2 = derivative(time + halfStepTime,nbSteps + halfH /*+halfH*k1*/,f_n,currentCelestialBody.velocity, currentCelestialBody,system ); // k2 computation
	for (int i=0; i < nbSteps; i++)
		f_n(i) = currentCelestialBody.position(i) + halfH*k2(i);
	k3 = derivative(time + halfStepTime,nbSteps + halfH /*+halfH*k2*/,f_n,currentCelestialBody.velocity,currentCelestialBody,system); // k3 computation
	for (int i=0; i < nbSteps; i++)
		f_n(i) = currentCelestialBody.position(i) + halfH*k3(i);
	k4 = derivative(time + halfStepTime,nbSteps + halfH /*+halfH*k3*/,f_n,currentCelestialBody.velocity,currentCelestialBody,system); // k4 computation
#pragma endregion
	// Then we have to set r. So yep, the new forces have to be calculated !!
	// Après chaque itération, ou après chhaque tour ? Chaque tour: sinon, on a des
	// Trucs relativement randoms!
	// And finally, we update the position vector/
	for (int i=0; i<nbSteps;i++)
		f_np1(i) = currentCelestialBody.position(i) + (double)(1/6)*h*(k1(i) + k2(i) + k3(i) + k4(i));
	currentCelestialBody.position = f_n;
	currentCelestialBody.velocity = g_n;
	// And lastly, we have to update our r in our global struct. Working with a pointer just to be sure that our new params are updated
	(*system).setForces();
	
	return;
}

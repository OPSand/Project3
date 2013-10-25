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
	cb->position = toCartesian2D(d, theta);
	cb->velocity = toCartesian2D(v, orthogonal2D(theta)); // counterclockwise
}
#pragma endregion

#pragma region Solver


vec derivative2D (const vec& posVel, CelestialBody* celestialBody, SolarSystem* system)
{
	vec derivatives = vec(4);// (0): x | (1): y | (2): v_x | (3): v_y
	vec dist = vec(2);
	dist.zeros();
	derivatives.zeros();
#pragma region Recalculating the forces
	int n = system->n();
	for (int i= 0; i< n; i++)
	{
		if (system->body(i)->name == (*celestialBody).name);
		else 
		{
			dist(0) = (system->body(i))->position(0) -posVel(0); // We re-evaluate the new r
			dist(1) = (system->body(i))->position(1) - posVel(1);
			vec r = (*celestialBody).position_diff(system->body(i)); // gives the force the proper direction ... 
			// We update the velocity
			derivatives(2) += (cG * (*celestialBody).mass * system->body(i)->mass / pow(dist(0), 3.0)) * r(0)*posVel(0); // Newton's law of gravity on x
			derivatives(3) += (cG * (*celestialBody).mass * system->body(i)->mass / pow(dist(1), 3.0)) * r(0)*posVel(1); // Newton's law of gravity on y
		}
	}
#pragma endregion
	//printf("%f || \t",derivatives(2));
	derivatives[0] = posVel(2); // Derivative of the position in x
	derivatives[1] = posVel(3); // Derivative of the position in y
	//printf("%f || \t", derivatives(0));
	return derivatives;
}

// So this is our Runge Kutta 4 algorithm in 2D for now. 
void rk4_2D(int dim, int h, int time,SolarSystem* system, CelestialBody *currentCelestialBody)
{
#pragma region About the first part
	//So to begin with, we'll only study a system with two points: the sun and the earth.
	// The method to do so is the same as the method to compute this with a lot more elements
	// The only thing changing is the flag set in the beginning of the project3.cpp main
#pragma endregion
	int nbEq = dim * 2; // We have to compute the velocity + position for all our dim.
	vec f_np1 = vec(nbEq);
	vec k1 = vec(nbEq);
	vec k2 = vec(nbEq);
	vec k3 = vec(nbEq);
	vec k4 = vec(nbEq);
	double halfH = 0.5*h;// We are computing the RG4 with half steps.
	vec f_n = vec(nbEq);
	for (int i= 0 ; i< dim; i++)
	{
		f_n(i) = (*currentCelestialBody).position(i); // (0): position x | (1): position y
		f_n(dim - 1 + i) = (*currentCelestialBody).velocity(i); // (2): velocity x | (3): velocity y
	}
#pragma region ki Computations
	k1 = derivative2D(f_n,currentCelestialBody,system); // k1 
	for (int i=0; i < dim; i++)
	{
		f_n(i) = (*currentCelestialBody).position(i) + halfH*k1(i);
		f_n(dim - 1 + i) = (*currentCelestialBody).velocity(i) + halfH*k1(dim - 1 + i);
	}
	k2 = derivative2D(f_n,currentCelestialBody,system ); // k2
	for (int i=0; i < dim; i++)
	{
		f_n(i) = (*currentCelestialBody).position(i) + halfH*k2(i);
		f_n(dim - 1 + i) = (*currentCelestialBody).velocity(i) + halfH*k2(dim - 1 + i);
	}
	k3 = derivative2D(f_n,currentCelestialBody,system); // k3
	for (int i=0; i < dim; i++)
	{
		f_n(i) = (*currentCelestialBody).position(i) + halfH*k3(i);
		f_n(dim - 1 + i) = (*currentCelestialBody).velocity(i) + halfH*k3(dim - 1 + i);
	}
	k4 = derivative2D(f_n,currentCelestialBody,system); // k4
#pragma endregion

	// And finally, we update the position vector/
	for (int i=0; i<dim;i++)
	{
		f_np1(i) = (*currentCelestialBody).position(i) + (double)(1/6)*h*(k1(i) + k2(i) + k3(i) + k4(i));
		f_np1(dim - 1 + i) = (*currentCelestialBody).velocity(i) + double(1/6)*h*(k1(dim - 1 + i) + k2(dim - 1 + i) + k3(dim - 1 + i) + k4(dim - 1 + i));
	}
	double ini_posi = (*currentCelestialBody).position(0); // XXX : To be removed
	printf("%f",ini_posi); // XXX : To be removed
	
	for(int i=0; i< dim; i++)
	{
		(*currentCelestialBody).position(i) = f_np1(i); 
		(*currentCelestialBody).velocity(i) = f_np1(dim - 1 + i);
	}

	double fini_posi = (*currentCelestialBody).position(0); // XXX : To be removed
	printf("%f:", fini_posi);// XXX : To be removed

	//system.setForces();
	
	return;
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
	const double DELTA_T = 24 * 60 * 60; // time step length (s)
	const int PLOT_EVERY = 1; // plot every ...th step
	const int N_PLOT = (N_STEPS / PLOT_EVERY); // how many steps we actually plot

	// compiler flags
	const bool ADD_JUPITER = true;
	const bool ADD_ALL = true; // include the other 7 planets
	const bool DEBUG = false; // use for debugging only
	const bool USE_EULER = false; // if true, use Runge-Kutta

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
		sun->velocity = (-system.totalMomentum() / sun->mass); // v = p/m;
		//assert( norm(system.totalMomentum(), DIM) == 0.0 ); // check that total momentum is actually 0
	}
#pragma endregion

#pragma region Debugging code (Here be dragons!)
	if( DEBUG )
	{
		system.setForces();
		for( int i = 0; i < system.n(); i++ )
		{
			cout << system.body(i)->name << endl << endl << "Force:" << endl << system.body(i)->force << endl << "Acc:" << endl << system.body(i)->acc() << endl << "Pos:" << endl << system.body(i)->position << endl << "Vel:" << endl << system.body(i)->velocity << endl;
		}

		/* 
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
		cout << jupiter->plot; */

		for( int i = 0; i < system.n(); i++ )
		{
			cout << i << ": " << system.body(i)->name << endl; // output planet names
		}
	}
#pragma endregion

#pragma region Iterate

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
				if( USE_EULER )
				{
					// acc -> velocity (Euler-Cromer, for testing only)
					cb->velocity += DELTA_T * cb->acc();

					// velocity -> position (Euler-Cromer, for testing only)
					cb->position += DELTA_T * cb->velocity;
				}
				else // use Runge-Kutta
				{
					// do stuffs
				}
			}

			if( i % PLOT_EVERY == 0 ) // we want to plot this step
			{
				cb->plotCurrentPosition();
			}
		}
	}
#pragma endregion

#pragma region Plot
	/* plot to file for each CB
	for( int j = 0; j < system.n(); j++ )
	{
		system.body(j)->positionToFile(); // saved as "<name>.dat"
	} */

	// plot X and Y coordinates for the entire system as matrices
	system.plotDim(0, "X.dat");
	system.plotDim(1, "Y.dat");

	cout << "Finished plotting " << N_PLOT << " of " << N_STEPS << " steps!";
#pragma endregion
	
#pragma region More debugging
	if ( DEBUG )
	{
		CelestialBody* Earth = system.body(1);
		int n = system.n();
		int t= 0;
		while (t < N_STEPS)
		{
			rk4_2D(2,DELTA_T,N_STEPS,&system,Earth);
			printf("t: % d | x: %f | y %f \t vx: %f | vy: %f",t,Earth->position(0),Earth->position(1),Earth->velocity(0),Earth->velocity(1));
			t+=DELTA_T;
		}
	}
#pragma endregion

	getchar(); // pause

	return 0; // exit
}


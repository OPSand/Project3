#include "stdafx.h"
#include "SolarSystem.h"

SolarSystem::SolarSystem(int dim)
{
	this->_dim = dim;
	this->_bodies = vector<CelestialBody>();
}


SolarSystem::~SolarSystem(void)
{
	// no uses of new
}


void SolarSystem::setForces(void)
{
	// gravitational constant
	const double G = 6.67385e-11; // N (m/kg)^2

	int n = this->n();
	for( int i = 0; i < (n - 1); i++ ) // i: 0 -> n-2
	{
		CelestialBody cb_i = this->body(i);
		
		// reset on first pass (0)
		if( i == 0 )
		{
			cb_i.force.fill(0);
		}

		for( int j = (i + 1); j < n; j++ ) // j: i+1 -> n-1
		{
			CelestialBody cb_j = this->body(j);

			// reset on first pass (1 -> n-1)
			if( i == 0 )
			{
				cb_j.force.fill(0);
			}
			
			double dist = cb_i.dist(cb_j); // distance (absolute value)
			vec r = cb_i.position_diff(cb_j); // gives the force the proper direction
			vec F = (G * cb_i.mass() * cb_j.mass() / pow(dist, 3.0)) * r; // Newton's law of gravity

			cb_i.force += F; // add force contribution to i
			cb_j.force -= F; // add force contribution to j (Newton's 3rd law)
		}
	}
}

// return a celestial body at the index i
CelestialBody SolarSystem::body(int i)
{
	return _bodies.at(i);
}

// add a new celecstial body to solar system
void SolarSystem::add(CelestialBody cb)
{
	this->_bodies.push_back(cb); // append to end of vector
}
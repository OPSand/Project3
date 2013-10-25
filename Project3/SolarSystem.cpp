#include "stdafx.h"
#include "SolarSystem.h"

SolarSystem::SolarSystem(int dim, int nSteps)
{
	this->_dim = dim;
	this->_nSteps = nSteps;
	this->_bodies = new vector<CelestialBody*>();
}


SolarSystem::~SolarSystem(void)
{
	while( this->n() > 0 ) // counter will be updated automatically
	{
		CelestialBody* cb = this->body(this->n() - 1); // return last element
		this->_bodies->pop_back(); // remove it from vector
		delete cb; // delete it
	}
	delete _bodies; // delete vector
}


void SolarSystem::setForces(void)
{
	int n = this->n();
	for( int i = 0; i < (n - 1); i++ ) // i: 0 -> n-2
	{
		CelestialBody* cb_i = this->body(i);
		
		// reset on first pass (0)
		if( i == 0 )
		{
			cb_i->force.fill(0);
		}

		for( int j = (i + 1); j < n; j++ ) // j: i+1 -> n-1
		{
			CelestialBody* cb_j = this->body(j);

			// reset on first pass (1 -> n-1)
			if( i == 0 )
			{
				cb_j->force.fill(0.0);
			}
			
			double dist = cb_i->dist(cb_j); // distance (absolute value)
			vec r = cb_i->position_diff(cb_j); // gives the force the proper direction
			vec F = (cG * cb_i->mass * cb_j->mass / pow(dist, 3.0)) * r; // Newton's law of gravity

			cb_i->force += F; // add force contribution to i
			cb_j->force -= F; // add force contribution to j (Newton's 3rd law)
		}
	}
}

// return a celestial body at the index i
CelestialBody* SolarSystem::body(int i)
{
	return _bodies->at(i);
}

// add a new celestial body to solar system
void SolarSystem::add(CelestialBody* cb)
{
	this->_bodies->push_back(cb); // append to end of vector
}

// return total momentum of system
vec SolarSystem::totalMomentum()
{
	vec mom(this->_dim);
	mom.fill(0.0);

	for(int i = 0; i < this->n(); i++)
	{
		CelestialBody* cb_i = this->body(i);
		if( ! cb_i->fixed ) // fixed bodies never move
		{
			mom += (cb_i->mass * cb_i->velocity); // add p = m*v (non-relativistic)
		}
	}

	return mom;
}
#include "stdafx.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

using namespace std;
using namespace arma;

// constructor
CelestialBody::CelestialBody(string name, double mass, SolarSystem system)
{
	// set protected member variables
	this->_name = name;
	this->_mass = mass;
	this->_dim = system.dim();

	// set dimensions
	this->position = vec(_dim);
	this->velocity = vec(_dim);
	this->force = vec(_dim);

	// default values
	this->fixed = false;
	this->position.fill(0);
	this->velocity.fill(0);
	this->force.fill(0);
}

// copy constructor
CelestialBody::CelestialBody(const CelestialBody &cb)
{
	this->_name = cb._name;
	this->_mass = cb._mass;
	this->_dim = cb._dim;

	this->fixed = cb.fixed;

	// Armadillo being Armadillo, these should be copied, not just passed by reference.
	this->position = cb.position;
	this->velocity = cb.velocity;
	this->force = cb.force;
}

// destructor
CelestialBody::~CelestialBody(void)
{
	// empty because we don't use new
}

// operator =
CelestialBody CelestialBody::operator = (const CelestialBody &cb)
{
	if( this != &cb ) // protect against invalid self-assignment
	{
		this->_name = cb._name;
		this->_mass = cb._mass;
		this->_dim = cb._dim;

		this->fixed = cb.fixed;

		// Armadillo being Armadillo, these should be copied, not just passed by reference.
		this->position = cb.position;
		this->velocity = cb.velocity;
		this->force = cb.force;
	}

	return *this; // convention
}
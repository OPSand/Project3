#include "stdafx.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

using namespace std;
using namespace arma;

// constructor (system passed by reference so we do not copy it)
CelestialBody::CelestialBody(const string& name, double mass, SolarSystem* system, bool fixed)
{
	// set protected member variables
	this->name = name;
	this->mass = mass;
	this->_dim = system->dim();

	// initialize vectors with correct dimension
	this->position = vec(this->_dim);
	this->velocity = vec(this->_dim);
	this->force = vec(this->_dim);

	// default values
	this->fixed = fixed;
	this->position.fill(0);
	this->velocity.fill(0);
	this->force.fill(0);

	// add to system
	system->add(this);
}

// copy constructor
CelestialBody::CelestialBody(const CelestialBody &cb)
{
	this->name = cb.name;
	this->mass = cb.mass;
	this->_dim = cb._dim;

	this->fixed = cb.fixed;

	// Armadillo being Armadillo, these should be copied, not just passed by reference (we hope).
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
		this->name = cb.name;
		this->mass = cb.mass;
		this->_dim = cb._dim;
		this->fixed = cb.fixed;

		// Armadillo being Armadillo, these should be copied, not just passed by reference (we hope).
		this->position = cb.position;
		this->velocity = cb.velocity;
		this->force = cb.force;
	}

	return *this; // to allow operator chaining: a = b = c
}
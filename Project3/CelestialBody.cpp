#include "stdafx.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

using namespace arma;

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
CelestialBody::CelestialBody(CelestialBody &cb)
{
	this->_name = cb._name;
	this->_mass = cb._mass;
	this->_dim = cb._dim;

	this->position = cb.position;
	//...
}


CelestialBody::~CelestialBody(void)
{
}

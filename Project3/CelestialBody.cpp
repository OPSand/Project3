#include "stdafx.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

CelestialBody::CelestialBody(string name, double mass, SolarSystem system)
{
	// set protected member variables
	_name = name;
	_mass = mass;
	_dim = system.dim();

	// set dimensions
	position(_dim);
	velocity(_dim);
	force(_dim);

	// default values
	fixed = false;
	position.fill(0);
	velocity.fill(0);
	force.fill(0);
}


CelestialBody::~CelestialBody(void)
{
}

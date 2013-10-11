#include "stdafx.h"
#include "CelestialBody.h"


CelestialBody::CelestialBody(string name, int dim, double mass)
{
	// set protected member variables
	_name = name;
	_mass = mass;
	_dim = dim;

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

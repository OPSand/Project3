#pragma once

#include "stdafx.h"
#include "SolarSystem.h"

using namespace arma;
using namespace std;

class SolarSystem; // forward declaration to avoid circular reference

class CelestialBody
{
public:
	CelestialBody(string name, double mass, SolarSystem system);
	CelestialBody(const CelestialBody &cb);
	~CelestialBody(void);
	CelestialBody operator = (const CelestialBody &cb);
	vec position;
	vec velocity;
	vec force;
	bool fixed;

protected:
	string _name;
	double _mass;
	int _dim;
public:

	string name(void)
	{
		return _name;
	}

	double mass(void)
	{
		return _mass;
	}

	// returns the acceleration when the force is set
	vec acc(void)
	{
		return (force/mass());
	}

	// returns the position of cb relative to this in vector form
	vec position_diff(CelestialBody cb)
	{
		return (cb.position - this->position);
	}

	// returns the distance between this and cb as a scalar
	double dist(CelestialBody cb)
	{
		return norm(this->position_diff(cb), _dim);
	}	
};


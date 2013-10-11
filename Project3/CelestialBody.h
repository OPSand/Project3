#pragma once

#include "stdafx.h"
#include "SolarSystem.h"

using namespace arma;
using namespace std;

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
public:

	string name(void)
	{
		return _name;
	}

	double mass(void)
	{
		return _mass;
	}

	vec acc(void)
	{
		return (force/mass());
	}

	double dist(CelestialBody cb)
	{
		return norm((position - cb.position), _dim);
	}
protected:
	int _dim;
};


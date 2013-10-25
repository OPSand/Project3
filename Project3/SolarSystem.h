#pragma once

#include "stdafx.h"
#include "CelestialBody.h"

class CelestialBody; // forward declaration to avoid circular reference

class SolarSystem
{
protected:
	int _dim;
	int _nSteps;
	vector<CelestialBody*>* _bodies; // list of celestial bodies in solar system (use pointers to avoid needless copying)

public:
	SolarSystem(int dim, int nSteps);
	~SolarSystem(void);
	void setForces(void);

	// return dimension of system
	int dim(void)
	{
		return _dim;
	}

	// return number of time steos
	int nSteps(void)
	{
		return _nSteps;
	}

	// return number of celestial bodies in system
	int n(void)
	{
		return _bodies->size();
	}

	// return a celestial body at the index i
	CelestialBody* body(int i);

	// add a new celecstial body to solar system
	void add(CelestialBody* cb);

	// return total momentum of system
	vec totalMomentum(void);
};
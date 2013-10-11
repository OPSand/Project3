#pragma once

#include "stdafx.h"

class SolarSystem
{
public:
	SolarSystem(void);
	~SolarSystem(void);
	void setForces(void);
//	list<CelestialBody> bodies;

	void add(CelestialBody cb)
	{
		// bodies.push_back(cb);
	}
protected:
	CelestialBody* _bodies;
//	int maxSize;
	int _capacity;
	int _count;
};


#pragma once

#include "CelestialBody.h"

class SolarSystem
{
public:
	SolarSystem(void);
	~SolarSystem(void);
	void setForces(void);
	list<CelestialBody> bodies;
};


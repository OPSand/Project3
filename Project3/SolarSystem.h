#pragma once

#include "stdafx.h"
#include <vector>

class SolarSystem
{
public:
	SolarSystem(int dim);
	~SolarSystem(void);
	void setForces(void);

protected:
	int _dim;
public:

	int dim(void)
	{
		return _dim;
	}
};
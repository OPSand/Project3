#include "stdafx.h"
#include "SolarSystem.h"


SolarSystem::SolarSystem(int dim)
{
	_dim = dim;

	// TODO: init bodies
}


SolarSystem::~SolarSystem(void)
{
}


void SolarSystem::setForces(void)
{
	// for each CB i
		// for each CB j > i
			// add force contribution F_ij and F_ji (=-F_ij) to EACH 
}

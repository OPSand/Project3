#pragma once
class CelestialBody
{
public:
	CelestialBody(string name, double mass);
	~CelestialBody(void);
	vec position;
	vec velocity;
	vec force;
	bool fixed;

	double acc(void)
	{
		return 0;
	}
protected:
	string _name;
	double _mass;
public:

	string name(void)
	{
		return string();
	}

	double mass(void)
	{
		return 0;
	}
};


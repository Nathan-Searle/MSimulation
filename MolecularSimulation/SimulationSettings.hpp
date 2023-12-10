#ifndef SIMULATIONSETTINGS_HPP
#define SIMULATIONSETTINGS_HPP

class SimulationSettings
{
public:
	// This is the major sin, but works for now (hardcoded)
	const static int totalIterations = 500;
	const static int nParticles = 21;
	constexpr static double LJEpsilon = 125.7 * 1.380649E-16 * 6.022E39; // amu * angstrom / (100 fs)^2
	constexpr static double LJSigma = 0.3345 * 10; // Angstrom
};

#endif
#ifndef SIMULATIONSETTINGS_HPP
#define SIMULATIONSETTINGS_HPP

class SimulationSettings
{
public:
	// This is the major sin, but works for now (hardcoded)
	const static int totalIterations = 150;
	const static int nParticles = 25;
	constexpr static double LJEpsilon = 125.7 * 1.380649E-16; // ergs
	constexpr static double LJSigma = 0.3345 * 10; // Angstrom
};

#endif
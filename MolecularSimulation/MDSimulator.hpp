#ifndef MDSIMULATOR_HPP
#define MDSIMULATOR_HPP
#include <array>
#include "SimulationSettings.hpp"

class MDSimulator
{
public:
	MDSimulator(double, double, double, double);
	void solve();

private:
	int currentIteration;
	int totalIteration;
	double length;
	double temperature;
	double timeDelta;
	double rcut;
	double mass;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100> positionArr;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100> velocityArr;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100> forceArr;

	void initialize();
	void velocityVerlet();
	const double potentialCalculator(double);
	const double kineticCalculator(double);
	double temperatureCalculator();
	void thermostat();
	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> findDistances();
};

#endif
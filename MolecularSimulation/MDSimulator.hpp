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
	double currentTime;
	double length;
	double temperature;
	double maxTime;
	double timeDelta;
	double rcut;
	double mass;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations> positionArr;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations> velocityArr;
	std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations> forceArr;

	void initialize();
	void velocityVerlet();
	double potentialCalculator(double);
	double kineticCalculator(double);
	double temperatureCalculator();
	void thermostat();
	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> findDistances();
	double forceCalculator(double);
};

#endif
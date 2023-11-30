#ifndef MDSIMULATOR_HPP
#define MDSIMULATOR_HPP
#include <array>

class MDSimulator
{
public:
	MDSimulator(double, double, double, double);
	void solve();

private:
	// This is the major sin, but works for now (hardcoded)
	const static int totalIterations = 100;
	const static int nParticles = 10;

	double currentTime;
	double length;
	double temperature;
	double maxTime;
	double timeDelta;
	double mass;
	std::array<std::array<std::array<double, 3>, nParticles>, totalIterations> positionArr;
	std::array<std::array<std::array<double, 3>, nParticles>, totalIterations> velocityArr;

	void velocityVerlet();
	double potentialCalculator(double);
	double kineticCalculator(double);
	// Need a thermometer here?
	void updatePositionArr();
};

#endif
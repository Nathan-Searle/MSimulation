#include "MDSimulator.hpp"

MDSimulator::MDSimulator(double particleMass, double lengthOfCube, double temperatureIn, double maxTime) :
	mass(particleMass),
	length(lengthOfCube),
	temperature(temperatureIn),
	currentTime(0),
	timeDelta(maxTime/totalIterations),
	positionArr(std::array<std::array<std::array<double, 3>, nParticles>, totalIterations>()),
	velocityArr(std::array<std::array<std::array<double, 3>, nParticles>, totalIterations>())
{
}

void MDSimulator::solve() {
	while (MDSimulator::currentTime < MDSimulator::totalIterations) {
		velocityVerlet();

		MDSimulator::currentTime++;
	}
}

void MDSimulator::velocityVerlet() {
	// implement 3D velocity verlet algorithm here
}

// Helpers of velocity verlet
double MDSimulator::potentialCalculator(double coordinate) {
	// implement the LJ potential here to return the value based on coordinate
	return 0.0;
}

double MDSimulator::kineticCalculator(double velocity) {
	// calculate the kinetic energy
	return 0.0;
}

void MDSimulator::updatePositionArr() {
	// Implement the periodic boundary condition
}
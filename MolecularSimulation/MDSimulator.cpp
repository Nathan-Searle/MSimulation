#include "MDSimulator.hpp"

MDSimulator::MDSimulator(int numberOfParticles, double particleMass, double lengthOfCube, double temperatureIn, double totalTime, int iterationCount) :
	positionArr(std::vector<double[3]>(numberOfParticles)),
	velocityArr(std::vector<double[3]>(numberOfParticles)),
	mass(particleMass),
	length(lengthOfCube),
	temperature(temperatureIn),
	maxTime(totalTime),
	totalIterations(iterationCount),
	currentTime(0),
	timeDelta(maxTime/totalIterations)
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
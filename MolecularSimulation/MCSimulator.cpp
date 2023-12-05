#include "MCSimulator.hpp"
#include "SimulationSettings.hpp"

MCSimulator::MCSimulator(double particleMass, double lengthOfCube, int maxIterationCount) :
	mass(particleMass),
	length(lengthOfCube),
	maxIterations(maxIterationCount)
{
}

void MCSimulator::solve() {
	initialize();

	runMetropolis();

	// Print out a bunch of stuff at the end here for properties, etc.
}

void MCSimulator::initialize() {
	double maxDistance = pow((pow(length, 3) / SimulationSettings::nParticles), 1.0 / 3.0);
	int nIntervals = (int)(length / maxDistance) + 1;
	int particleIndex;
	for (int i = 0; i < nIntervals; i++) {
		for (int j = 0; j < nIntervals; j++) {
			for (int k = 0; k < nIntervals; k++) {
				particleIndex = (int)(j * nIntervals + k + i * nIntervals * nIntervals);
				if (particleIndex == SimulationSettings::nParticles)
					return;
				//positionArr[0][particleIndex][0] = -length / 2 + maxDistance / 2 + maxDistance * i;
				//positionArr[0][particleIndex][1] = -length / 2 + maxDistance / 2 + maxDistance * j;
				//positionArr[0][particleIndex][2] = -length / 2 + maxDistance / 2 + maxDistance * k;
			}
		}
	}
}


void MCSimulator::runMetropolis() {

}
#include "MDSimulator.hpp"
#include <iostream>
#include <fstream>

MDSimulator::MDSimulator(double particleMass, double lengthOfCube, double temperatureIn, double maxTime) :
	mass(particleMass),
	length(lengthOfCube),
	temperature(temperatureIn),
	rcut(lengthOfCube / 2),
	currentIteration(0),
	totalIteration(0),
	timeDelta(maxTime / SimulationSettings::totalIterations),
	positionArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100>()),
	velocityArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100>()),
	forceArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, 100>())
{
}

void MDSimulator::solve() {
	initialize();

	std::ofstream tradjectories;
	tradjectories.open("trajectories.txt");

	while (totalIteration < SimulationSettings::totalIterations - 1) {
		velocityVerlet();

		if (currentIteration % 99 == 98) {
			// Save some results periodically
			for (int j = 0; j < currentIteration+1; j++) {
				tradjectories << "[";
				for (int i = 0; i < SimulationSettings::nParticles; i++) {
					tradjectories << "[" << positionArr[j][i][0] << "," << positionArr[j][i][1] << "," << positionArr[j][i][2] << "]";
					if (i <= SimulationSettings::nParticles - 2)
						tradjectories << ",";
				}
				if (j < positionArr.size())
					tradjectories << "],";
			}

			// Reset where things are (save memory)
			for (int i = 0; i < SimulationSettings::nParticles; i++) {
				for (int j = 0; j < 3; j++) {
					positionArr[0][i][j] = positionArr[98][i][j];
					positionArr[1][i][j] = positionArr[99][i][j];

					velocityArr[0][i][j] = velocityArr[98][i][j];
					velocityArr[1][i][j] = velocityArr[99][i][j];

					forceArr[0][i][j] = forceArr[98][i][j];
					forceArr[1][i][j] = forceArr[99][i][j];
				}
			}

			currentIteration = -1;
		}

		currentIteration++;
		totalIteration++;

		// thermostat();
	}

	// Save the remaining results
	for (int j = 0; j < totalIteration % 99 + 1; j++) {
		tradjectories << "[";
		for (int i = 0; i < SimulationSettings::nParticles; i++) {
			tradjectories << "[" << positionArr[j][i][0] << "," << positionArr[j][i][1] << "," << positionArr[j][i][2] << "]";
			if (i <= SimulationSettings::nParticles - 2)
				tradjectories << ",";
		}
		if (j < positionArr.size() - 1)
			tradjectories << "],";
		else
			tradjectories << "]";
	}
	
	tradjectories.close();
}

void MDSimulator::initialize() {
	double maxDistance = pow((pow(length, 3) / SimulationSettings::nParticles), 1.0 / 3.0);
	int nIntervals = (int)(length / (maxDistance / 1.1));
	double initVelocity = sqrt(temperature * 3 * 1.380649E-16 * 6.022E39 / mass); // angstrom / [time]
	int particleIndex;
	for (int i = 0; i < nIntervals; i++) {
		for (int j = 0; j < nIntervals; j++) {
			for (int k = 0; k < nIntervals; k++) {
				particleIndex = (int)(j * nIntervals + k + i * nIntervals * nIntervals);
				if (particleIndex == SimulationSettings::nParticles)
					goto calculateInitialForces;
				positionArr[0][particleIndex][0] = -length / 2 + maxDistance / 2 + maxDistance * i;
				positionArr[0][particleIndex][1] = -length / 2 + maxDistance / 2 + maxDistance * j;
				positionArr[0][particleIndex][2] = -length / 2 + maxDistance / 2 + maxDistance * k;

				int coordinate = rand() % 100;
				if (coordinate < 33)
					velocityArr[0][particleIndex][0] = initVelocity; // x

				else if (coordinate < 66)
					velocityArr[0][particleIndex][1] = initVelocity; // y

				else
					velocityArr[0][particleIndex][2] = initVelocity; // z

				if (rand() % 100 < 50) { // make negative "randomly"
					velocityArr[0][particleIndex][0] *= -1;
					velocityArr[0][particleIndex][1] *= -1;
					velocityArr[0][particleIndex][2] *= -1;
				}
			}
		}
	}

	calculateInitialForces:
	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> allDistances = findDistances();
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// Start at 0 force in each dimension
		forceArr[0][i][0] = 0;
		forceArr[0][i][1] = 0;
		forceArr[0][i][2] = 0;
		std::array < std::array<double, 3>, SimulationSettings::nParticles> distanceArr = allDistances[i];
		for (int j = 0; j < size(allDistances); j++) {
			// Find the total distance between each particle
			double r = sqrt(pow(distanceArr[j][0], 2) + pow(distanceArr[j][1], 2) + pow(distanceArr[j][2], 2));
			if (r > rcut) {
				forceArr[0][i][0] += 0;
				forceArr[0][i][1] += 0;
				forceArr[0][i][2] += 0;
			}
			else {
				double force = -24 * SimulationSettings::LJEpsilon * (pow(SimulationSettings::LJSigma, 6) / pow(r, 7) - 2 * pow(SimulationSettings::LJSigma, 12) / pow(r, 13));
				forceArr[0][i][0] += force * distanceArr[j][0];
				forceArr[0][i][1] += force * distanceArr[j][1];
				forceArr[0][i][2] += force * distanceArr[j][2];
			}
		}
	}
}

void MDSimulator::velocityVerlet() {
	// Update the positions
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// For each dimension
		for (int j = 0; j < 3; j++) {
			double nextPos = positionArr[currentIteration][i][j] + velocityArr[currentIteration][i][j] * timeDelta + pow(timeDelta, 2) / 2 * (velocityArr[currentIteration][i][j] / mass);

			// Break if something has gone horribly wrong
			if (abs(nextPos) > length) {
				std::cout << "We done lost a particle." << std::endl;
			}

			// Enforce periodic boundary condition
			if (nextPos > length / 2)
				nextPos -= length;
			else if (nextPos < -length / 2)
				nextPos += length;
			positionArr[currentIteration + 1][i][j] = nextPos;
		}
	}

	// Update the forces
	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> allDistances = findDistances();
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// Start at 0 force in each dimension
		forceArr[currentIteration + 1][i][0] = 0;
		forceArr[currentIteration + 1][i][1] = 0;
		forceArr[currentIteration + 1][i][2] = 0;
		std::array < std::array<double, 3>, SimulationSettings::nParticles> distanceArr = allDistances[i];
		for (int j = 0; j < allDistances.size(); j++) {
			// Find the total distance between each particle
			double r = sqrt(pow(distanceArr[j][0], 2) + pow(distanceArr[j][1], 2) + pow(distanceArr[j][2], 2));
			if (r > rcut) {
				forceArr[currentIteration + 1][i][0] += 0;
				forceArr[currentIteration + 1][i][1] += 0;
				forceArr[currentIteration + 1][i][2] += 0;
			}
			else {
				double force = -24 * SimulationSettings::LJEpsilon * (pow(SimulationSettings::LJSigma, 6) / pow(r, 7) - 2 * pow(SimulationSettings::LJSigma, 12) / pow(r, 13));
				forceArr[currentIteration + 1][i][0] += force * distanceArr[j][0];
				forceArr[currentIteration + 1][i][1] += force * distanceArr[j][1];
				forceArr[currentIteration + 1][i][2] += force * distanceArr[j][2];
			}
		}
	}

	// Update the velocities
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// For each dimension
		for (int j = 0; j < 3; j++) {
			velocityArr[currentIteration + 1][i][j] = velocityArr[currentIteration][i][j] + timeDelta / (2 * mass) * (forceArr[currentIteration][i][j] + forceArr[currentIteration + 1][i][j]);
		}
	}
}

const double MDSimulator::potentialCalculator(double coordinate) {
	double potential;
	if (coordinate > rcut)
		potential = 0;
	else {
		double sigmaGoo = pow(SimulationSettings::LJSigma / coordinate, 12) - pow(SimulationSettings::LJSigma / coordinate, 6);
		potential = 4 * SimulationSettings::LJEpsilon * sigmaGoo;
	}
	return potential; // ergs
}

const double MDSimulator::kineticCalculator(double velocity) {
	return 1.0 / 2.0 * mass * pow(velocity, 2); // ergs
}

double MDSimulator::temperatureCalculator() {
	//double initVelocity = sqrt(temperature * 3 * 1.380649E-16 * 6.022E39 / mass);
	const double firstTerm = 2 / (3 * 1.380649E-16 * 6.022E39 * SimulationSettings::nParticles);
	double sum = 0;
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = 0; j < 3; j++) {
			sum += kineticCalculator(velocityArr[currentIteration][i][j]);
		}
	}
	return firstTerm * sum; // K
}

std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> MDSimulator::findDistances() {

	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> allDistances;
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = i + 1; j < SimulationSettings::nParticles; j++) {
			double xDistance = positionArr[currentIteration][i][0] - positionArr[currentIteration][j][0]; // xi - xj
			double yDistance = positionArr[currentIteration][i][1] - positionArr[currentIteration][j][1]; // yi - yj
			double zDistance = positionArr[currentIteration][i][2] - positionArr[currentIteration][j][2]; // zi - zj

			// Check periodic boundary conditions for distance calculations
			if (xDistance < (-length / 2))
				xDistance += length;
			else if (xDistance > (length / 2))
				xDistance -= length;

			if (yDistance < (-length / 2))
				yDistance += length;
			else if (yDistance > (length / 2))
				yDistance -= length;

			if (zDistance < (-length / 2))
				zDistance += length;
			else if (zDistance > (length / 2))
				zDistance -= length;

			allDistances[i][j][0] = xDistance;
			allDistances[i][j][1] = yDistance;
			allDistances[i][j][2] = zDistance;

			allDistances[j][i][0] = -xDistance;
			allDistances[j][i][1] = -yDistance;
			allDistances[j][i][2] = -zDistance;
		}
	}
	return allDistances;
}

void MDSimulator::thermostat() {
	double currentTemp = temperatureCalculator();
	std::cout << currentTemp << ", ";
	double correction = 1.022E23 * (currentTemp - temperature);
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = 0; j < 3; j++) {
			forceArr[currentIteration][i][j] -= correction;
		}
	}
}
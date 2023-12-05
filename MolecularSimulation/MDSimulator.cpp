#include "MDSimulator.hpp"

MDSimulator::MDSimulator(double particleMass, double lengthOfCube, double temperatureIn, double maxTime) :
	mass(particleMass),
	length(lengthOfCube),
	temperature(temperatureIn),
	rcut(lengthOfCube / 2 / 1.5),
	currentTime(0),
	timeDelta(maxTime / SimulationSettings::totalIterations),
	positionArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations>()),
	velocityArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations>()),
	forceArr(std::array<std::array<std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::totalIterations>())
{
}

void MDSimulator::solve() {
	initialize();

	while (currentTime < SimulationSettings::totalIterations - 1) {
		velocityVerlet();

		thermostat();

		currentTime++;
	}

	// Print out a bunch of stuff at the end here for properties, etc.
}

void MDSimulator::initialize() {
	double maxDistance = pow((pow(length, 3) / SimulationSettings::nParticles), 1.0 / 3.0);
	int nIntervals = (int)(length / maxDistance) + 1;
	double initVelocity = 0;
	int particleIndex;
	for (int i = 0; i < nIntervals; i++) {
		for (int j = 0; j < nIntervals; j++) {
			for (int k = 0; k < nIntervals; k++) {
				particleIndex = (int)(j * nIntervals + k + i * nIntervals * nIntervals);
				if (particleIndex == SimulationSettings::nParticles)
					return;
				positionArr[0][particleIndex][0] = -length / 2 + maxDistance / 2 + maxDistance * i;
				positionArr[0][particleIndex][1] = -length / 2 + maxDistance / 2 + maxDistance * j;
				positionArr[0][particleIndex][2] = -length / 2 + maxDistance / 2 + maxDistance * k;

				int coordinate = rand() % 100;
				if (coordinate < 33) // x
					velocityArr[0][particleIndex][0] = initVelocity;

				else if (coordinate < 66) // y
					velocityArr[0][particleIndex][1] = initVelocity;

				else // z
					velocityArr[0][particleIndex][2] = initVelocity;

				if (rand() % 100 < 50) // make negative
					velocityArr[0][particleIndex][0] *= -1;
				velocityArr[0][particleIndex][1] *= -1;
				velocityArr[0][particleIndex][2] *= -1;
			}
		}
	}
}

void MDSimulator::velocityVerlet() {
	// Update the positions
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// For each dimension
		for (int j = 0; j < 3; j++) {
			double nextPos = positionArr[currentTime][i][j] + velocityArr[currentTime][i][j] * timeDelta + pow(timeDelta, 2) / 2 * (velocityArr[currentTime][i][j] / mass);

			// Enforce periodic boundary condition
			if (nextPos > length / 2)
				nextPos -= length;
			else if (nextPos < -length / 2)
				nextPos += length;
			positionArr[currentTime + 1][i][j] = nextPos;
		}
	}

	// Update the forces
	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> allDistances = findDistances();
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// Start at 0 force in each dimension
		forceArr[currentTime + 1][i][0] = 0;
		forceArr[currentTime + 1][i][1] = 0;
		forceArr[currentTime + 1][i][2] = 0;
		std::array < std::array<double, 3>, SimulationSettings::nParticles> distanceArr = allDistances[i];
		for (int j = 0; j < size(allDistances); j++) {
			// Find the total distance between each particle
			double r = sqrt(pow(distanceArr[j][0], 2) + pow(distanceArr[j][1], 2) + pow(distanceArr[j][2], 2));
			if (r > rcut) {
				forceArr[currentTime + 1][i][0] += 0;
				forceArr[currentTime + 1][i][1] += 0;
				forceArr[currentTime + 1][i][2] += 0;
			}
			else {
				double force = forceCalculator(r);
				forceArr[currentTime + 1][i][0] += force * distanceArr[j][0];
				forceArr[currentTime + 1][i][0] += force * distanceArr[j][1];
				forceArr[currentTime + 1][i][0] += force * distanceArr[j][2];
			}
		}
	}

	// Update the velocities
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		// For each dimension
		for (int j = 0; j < 3; j++) {
			velocityArr[currentTime + 1][i][j] = velocityArr[currentTime][i][j] + timeDelta / (2 * mass) * (forceArr[currentTime][i][j] + forceArr[currentTime + 1][i][j]);
		}
	}
}

double MDSimulator::potentialCalculator(double coordinate) {
	double potential;
	if (coordinate > rcut)
		potential = 0;
	else {
		double sigmaGoo = pow(SimulationSettings::LJSigma / coordinate, 12) - pow(SimulationSettings::LJSigma / coordinate, 6);
		potential = 4 * SimulationSettings::LJEpsilon * sigmaGoo;
	}
	return potential; // ergs
}

double MDSimulator::kineticCalculator(double velocity) {
	return 1 / 2 * mass * pow(velocity, 2); // ergs
}

double MDSimulator::temperatureCalculator() {
	const double firstTerm = 2 / (3 * 1.380649E-16 * SimulationSettings::nParticles);
	double sum = 0;
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = 0; j < 3; j++) {
			sum += kineticCalculator(velocityArr[currentTime][i][j]);
		}
	}
	return sum; // K
}

double MDSimulator::forceCalculator(double r) {
	return -24 * SimulationSettings::LJEpsilon * (pow(SimulationSettings::LJSigma, 6) / pow(r, 7) - 2 * pow(SimulationSettings::LJSigma, 12) / pow(r, 13));
}

std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> MDSimulator::findDistances() {

	std::array<std::array < std::array<double, 3>, SimulationSettings::nParticles>, SimulationSettings::nParticles> allDistances;
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = i + 1; j < SimulationSettings::nParticles; j++) {
			double xDistance = positionArr[currentTime][i][0] - positionArr[currentTime][j][0]; // xi - xj
			double yDistance = positionArr[currentTime][i][1] - positionArr[currentTime][j][1]; // yi - yj
			double zDistance = positionArr[currentTime][i][2] - positionArr[currentTime][j][2]; // zi - zj

			if (xDistance < (-length / 2))
				xDistance += length;

			allDistances[i][j][0] = xDistance;
			allDistances[i][j][1] = yDistance;
			allDistances[i][j][2] = zDistance;

			allDistances[j][i][0] = xDistance;
			allDistances[j][i][1] = yDistance;
			allDistances[j][i][2] = zDistance;
		}
	}
	return allDistances;
}

void MDSimulator::thermostat() {
	double currentTemp = temperatureCalculator(); // Use this value?
	for (int i = 0; i < SimulationSettings::nParticles; i++) {
		for (int j = 0; j < 3; j++) {
			forceArr[currentTime][i][j] += 0; // TODO : Need to add dummy force here to keep temperature good
		}
	}
}
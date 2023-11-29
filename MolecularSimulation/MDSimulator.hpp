#ifndef MDSIMULATOR_HPP
#define MDSIMULATOR_HPP
#include <vector>

class MDSimulator
{
public:
	MDSimulator(int, double, double, double, double, int);
	void solve();

private:
	double currentTime;
	std::vector<double[3]> positionArr;
	std::vector<double[3]> velocityArr;
	int totalIterations;
	double length;
	double temperature;
	double maxTime;
	double timeDelta;
	double mass;

	void velocityVerlet();
	double potentialCalculator(double);
	double kineticCalculator(double);
	// Need a thermometer here?
	void updatePositionArr();
};

#endif
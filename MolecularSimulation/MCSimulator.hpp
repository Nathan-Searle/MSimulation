#ifndef MCSIMULATOR_HPP
#define MCSIMULATOR_HPP
#include <array>

class MCSimulator
{
public:
	MCSimulator(double, double, int);
	void solve();

private:
	double mass;
	double length;
	int maxIterations;

	void initialize();
	void runMetropolis();
};

#endif
#include "MDSimulator.hpp"
#include "MCSimulator.hpp"
#include "SimulationSettings.hpp"
#include <array>
#include <iostream>

int main()
{
    //// setup number of particles and total number of iterations first, so this thing knows how to initialize memory

    //std::array<std::array<std::array<double, 3>, 10>, 100> posTime = std::array<std::array<std::array<double, 3>, 10>, 100>();

    //int i = 0;

    MDSimulator dynamicSimulator = MDSimulator(39.948, 1, 50, SimulationSettings::totalIterations * 0.05); // amu, angstrom, K, angstrom * amu^-1/2 * erg^1/2
    dynamicSimulator.solve();

    //MCSimulator metroSimulator = MCSimulator(39.948, 10, 1000); // amu, angstrom, number (of iterations allowed)
    //metroSimulator.solve();
}
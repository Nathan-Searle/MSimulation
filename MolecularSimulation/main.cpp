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

    MDSimulator dynamicSimulator = MDSimulator(39.948, 10, 85, SimulationSettings::totalIterations * 1E-16); // amu, angstrom, K, 1 femptosecond per iteration
    dynamicSimulator.solve();

    //MCSimulator metroSimulator = MCSimulator(39.948, 10, 1000); // amu, angstrom, number (of iterations allowed)
    //metroSimulator.solve();
}
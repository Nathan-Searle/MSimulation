//#include "MDSimulator.hpp"
#include <array>
#include <iostream>

int main()
{
    // setup number of particles and total number of iterations first, so this thing knows how to initialize memory

    std::array<std::array<std::array<double, 3>, 10>, 100> posTime = std::array<std::array<std::array<double, 3>, 10>, 100>();

    int i = 0;

    //MDSimulator dynamicSimulator = MDSimulator(3, 15, 3.2, 273.15, 100, 1000);
    //dynamicSimulator.solve();
}
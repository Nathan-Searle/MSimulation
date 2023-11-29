#include <iostream>
#include "VelocityVerlet.hpp"

int main()
{
    std::cout << "Hello World!\n" << std::endl;
    VelocityVerlet solver = VelocityVerlet();
    std::cout << "Bana_Func (9) = " << solver.bana_func(9);
}
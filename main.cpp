#include "new_core.hpp"

int main(int ac, char **av)
{
    try {
        algo::Core core;
        core.main(ac, av);
    } catch(std::exception &e) {
        std::cout << e.what() << std::endl;
    }
    return (0);
}
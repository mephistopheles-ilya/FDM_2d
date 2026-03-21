#include <iostream>
#include "parse_command_line.hpp"


int main(int argc, char **argv)
{
    Parser parser;
    parser. template add<double> ("hx", 0.1);
    parser. template add<double> ("hy", 0.1);
    parser. template add<double> ("ht", 0.1);
    parser. template add<double> ("mu", 0.1);
    parser. template add<double> ("pp", 1.);
    parser. template add<unsigned int> ("Nx", 100);
    parser. template add<unsigned int> ("Ny", 100);
    parser. template add<unsigned int> ("Mt", 100);
    parser. template add<unsigned int> ("k", 0, "nested parametr 2^k");
    parser. template add<unsigned int> ("n", 1, "number of threads");
    parser. template add<bool> ("Laspack", false, "use Laspacke for solver");
    parser. template add<bool> ("Eigen", false, "use Eigen for solver");

    int ret = 0;
    ret = parser.parse(argc, argv);
    if (ret < 0)
      return 0;
   
    parser.print_values ();

    

    return 0;
}

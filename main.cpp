#include <iostream>
#include "matrix_storage.hpp"
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
    parser. template add<unsigned int> ("Nt", 100);
    parser. template add<unsigned int> ("k", 0, "nested parametr 2^k");
    parser. template add<unsigned int> ("n", 1, "number of threads");
    parser. template add<bool> ("Laspack", false, "use Laspacke for solver");
    parser. template add<bool> ("Eigen", false, "use Eigen for solver");

    int ret = 0;
    ret = parser.parse (argc, argv);
    if (ret < 0)
      return ret;
   
    parser.print_values ();

    matrix_storage matrix_rhs;
    ret = matrix_rhs.init_grid (parser);
    if (ret < 0)
      return ret;

    ret = matrix_rhs.set_diff_params (parser);
    if (ret < 0)
      return ret;

    ret = matrix_rhs.allocate ();
    if (ret < 0)
      return ret;

    ret = matrix_rhs.fill_matrix_pattern ();
    if (ret < 0)
      return ret;

    unsigned int Nt = 0;
    ret = parser.get ("Nt", Nt);
    if (ret < 0)
      return ret;

    unsigned int step = 0;
    for (step = 0; step <= Nt; ++step)
      {
        matrix_rhs.fill_matrix (step);
        matrix_rhs.solve ();
      }

    matrix_rhs.deallocate ();

    return 0;
}

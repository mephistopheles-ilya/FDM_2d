#include "matrix_storage.hpp"
#include "parse_command_line.hpp"
#include <fenv.h>



int main(int argc, char **argv)
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    Parser parser;
    parser. template add<double> ("hx", 0.01);
    parser. template add<double> ("hy", 0.01);
    parser. template add<double> ("ht", 0.01);

    parser. template add<double> ("mu", 0.1);
    parser. template add<double> ("pp", 1.);

    parser. template add<double> ("eps", 1e-8);
    parser. template add<unsigned int> ("maxit", 2000);

    parser. template add<unsigned int> ("solver", 0);

    parser. template add<unsigned int> ("np", 0);

    int ret = 0;
    ret = parser.parse (argc, argv);
    if (ret < 0)
      return ret;
   
    parser.print_values ();

    matrix_storage matrix_rhs;
    ret = matrix_rhs.prepare_computations (parser, 0);
    if (ret < 0)
      return ret;

    ret = matrix_rhs.compute_solution (nullptr, 0);
    if (ret < 0)
      return ret;

    unsigned int np = 0;
    parser.get ("np", np);
    for (unsigned int k = 1; k <= np; ++k)
      {
        matrix_storage nested_matrix_rhs;
        ret = nested_matrix_rhs.prepare_computations (parser, k);
        if (ret < 0)
          return ret;

        ret = nested_matrix_rhs.compute_solution (&matrix_rhs, k);
        if (ret < 0)
          return ret;
      }

    return 0;
}

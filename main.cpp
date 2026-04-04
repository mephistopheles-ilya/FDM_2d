#include <iostream>
#include "matrix_storage.hpp"
#include "parse_command_line.hpp"
#include <fenv.h>

int main(int argc, char **argv)
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    Parser parser;
    parser. template add<double> ("hx", 0.01);
    parser. template add<double> ("hy", 0.01);
    parser. template add<double> ("ht", 0.01);
    parser. template add<double> ("mu", 0.1);
    parser. template add<double> ("pp", 1.);
    parser. template add<unsigned int> ("Nx", 100);
    parser. template add<unsigned int> ("Ny", 100);
    parser. template add<unsigned int> ("Nt", 100);
    parser. template add<double> ("eps", 1e-3);
    parser. template add<int> ("maxit", 100);

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

    ret = matrix_rhs.init_solver (parser);
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
        ret = matrix_rhs.solve ();
        if (ret < 0)
          return ret;
        double C_norm_G = matrix_rhs. template calculate_C_norm <G> (step);
        double C_norm_V1 = matrix_rhs. template calculate_C_norm <V1> (step);
        double C_norm_V2 = matrix_rhs. template calculate_C_norm <V2> (step);

        double L2_norm_G = matrix_rhs. template calculate_L2_norm <G> (step);
        double L2_norm_V1 = matrix_rhs. template calculate_L2_norm <V1> (step);
        double L2_norm_V2 = matrix_rhs. template calculate_L2_norm <V2> (step);

        double W1_norm_G = matrix_rhs. template calculate_W1_norm <G> (step);
        double W1_norm_V1 = matrix_rhs. template calculate_W1_norm <V1> (step);
        double W1_norm_V2 = matrix_rhs. template calculate_W1_norm <V2> (step);

        printf ("Time step: %d\n", step);
        printf("C_nrom: G = %e, V1 = %e, V2 = %e\n", C_norm_G, C_norm_V1, C_norm_V2); 
        printf("L2_nrom: G = %e, V1 = %e, V2 = %e\n", L2_norm_G, L2_norm_V1, L2_norm_V2); 
        printf("W1_nrom: G = %e, V1 = %e, V2 = %e\n", W1_norm_G, W1_norm_V1, W1_norm_V2); 

        matrix_rhs.update_prev_solution ();
      }

    return 0;
}

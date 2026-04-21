#include "matrix_storage.hpp"

template <unsigned int solver_type>
unsigned int matrix_storage::fill_matrix (unsigned int time_step)
  {
    unsigned int n_elements = grid.get_n_elements ();
    double hx = 0;
    double hy = 0;
    double ht = 0;
    grid.get_h (&hx, &hy, &ht);

    double min_g = min_G ();
    double mum = mu * std::exp (-min_g);
    for (unsigned int element_i = 0; element_i < n_elements; element_i++)
      {
        unsigned int i = 0;
        unsigned int j = 0;
        grid.convert_element_i_to_ij (element_i, i, j);
        unsigned int border_type = grid.get_bored_type (i, j);

        switch (border_type)
        {
          // inner domain
          case INNER:
            {
              // equation for G
              set_diag <solver_type> (4, G, element_i);
              set_off_diag<solver_type> (-ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)), G, i - 1, j, element_i, G); 
              set_off_diag<solver_type> (-ht / hy * (GVVn (V2, i, j - 1) + GVVn (V2, i, j)), G, i, j - 1, element_i, G);
              set_off_diag<solver_type> (ht / hx  * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)), G, i + 1, j, element_i, G);
              set_off_diag<solver_type> (ht / hy  * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)), G, i, j + 1, element_i, G);

              set_off_diag<solver_type> (-2 * ht / hx, V1, i - 1, j, element_i, G);
              set_off_diag<solver_type> (-2 * ht / hy, V2, i, j - 1, element_i, G);
              set_off_diag<solver_type> (2 * ht / hx, V1, i + 1, j, element_i, G);
              set_off_diag<solver_type> (2 * ht / hy, V2, i, j + 1, element_i, G);

              set_rhs (4 * GVVn (G, i, j) + ht * GVVn (G, i, j) * 
                  (1. / hx * (GVVn (V1, i + 1, j) - GVVn (V1, i - 1, j)) + 1. / hy * (GVVn (V2, i, j + 1) - GVVn (V2, i, j - 1))) + 4 * ht * Func_0 (time_step * ht, i * hx, j * hy) , G, element_i);

              // equation for V1
              set_diag<solver_type> (6 + 4 * ht * mum * (4 / hx / hx + 3 / hy / hy), V1, element_i);
              set_off_diag<solver_type> (- (ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)) + 8 * ht * mum / hx / hx), V1, i - 1, j, element_i, V1);
              set_off_diag<solver_type> (-(3 * ht / 2 / hy * (GVVn (V2, i, j - 1) + GVVn (V2 , i, j)) + 6 * ht * mum / hy / hy), V1, i, j - 1, element_i, V1);
              set_off_diag<solver_type> (ht / hx * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)) - 8 * ht * mum / hx / hx, V1, i + 1, j, element_i, V1);
              set_off_diag<solver_type> (3 * ht / 2 / hy * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)) - 6 * ht * mum / hy / hy, V1, i, j + 1, element_i, V1);

              double H = std::exp (GVVn (G, i, j));
              set_off_diag<solver_type> (-3 * ht * pd (H, pp) / hx, G, i - 1, j, element_i, V1);
              set_off_diag<solver_type> (3 * ht * pd (H, pp) / hx, G, i + 1, j, element_i, V1);

              set_rhs (6 * GVVn (V1, i, j) + 3 * ht / 2 / hy * GVVn (V1, i, j) * (GVVn (V2, i, j + 1) - GVVn (V2, i, j - 1)) + 6 * ht * (mu / H - mum) * (4. / 3 / hx / hx * (GVVn (V1, i + 1, j)
                        - 2 * GVVn (V1, i, j) + GVVn (V1, i - 1, j)) + 1 / hy / hy * (GVVn ( V1, i, j + 1) - 2 * GVVn (V1, i, j) + GVVn (V1, i, j - 1))) + 
                    ht * mu / 2 / H / hx / hy * (GVVn (V2, i + 1, j + 1) - GVVn (V2, i - 1, j + 1) - GVVn (V2, i + 1, j - 1) + GVVn (V2, i - 1, j - 1)) + 
                    6 * ht * Func_1 (time_step * ht, i * hx, j * hy, pp, mu)
                    , V1, element_i);

              // equation for V2
              set_diag<solver_type> (6 + 4 * ht * mum * (3 / hx / hx + 4 / hy / hy), V2, element_i);
              set_off_diag<solver_type> (-(3 * ht / 2 / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)) + 6 * ht * mum / hx / hx), V2, i - 1, j, element_i, V2);
              set_off_diag<solver_type> (-(ht / hy * (GVVn (V2, i, j - 1) + GVVn (V2, i, j)) + 8 * ht * mum / hy / hy), V2, i, j - 1, element_i, V2);
              set_off_diag<solver_type> (3 * ht / 2 / hx * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)) - 6 * ht * mum / hx / hx, V2, i + 1, j, element_i, V2);
              set_off_diag<solver_type> (ht / hy * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)) - 8 * ht * mum / hy / hy, V2, i, j + 1, element_i, V2);

              set_off_diag<solver_type> (-3 * ht * pd (H, pp) / hy, G, i, j - 1, element_i, V2);
              set_off_diag<solver_type> (3 * ht * pd (H, pp) / hy, G, i, j + 1, element_i, V2);

              set_rhs (6 * GVVn (V2, i, j) + 3 * ht / 2 / hx * GVVn (V2, i, j) * (GVVn (V1, i + 1, j) - GVVn (V1, i - 1, j)) + 6 * ht * (mu / H - mum) * (1 / hx / hx *
                    (GVVn (V2, i + 1, j) - 2 * GVVn (V2, i, j) + GVVn (V2, i - 1, j)) + 4. / 3 / hy / hy * (GVVn (V2, i, j + 1) - 2 * GVVn (V2, i, j) + GVVn (V2, i, j - 1))) + 
                  ht * mu / 2 / H / hx / hy * (GVVn (V1, i + 1, j + 1) - GVVn (V1, i - 1, j + 1) - GVVn (V1, i + 1, j - 1) + GVVn (V1, i - 1, j - 1)) + 
                  6 * ht * Func_2 (time_step * ht, i * hx, j * hy, pp, mu)
                  , V2, element_i);

                break;
              }

            // left boundary
            case X_LEFT:
              {
                // equation for G
                set_diag<solver_type> (2 - ht/hx * GVVn (V1, i, j), G, element_i);
                set_off_diag <solver_type> (ht/hx * GVVn (V1, i + 1, j), G, i + 1, j, element_i, G);
                set_off_diag <solver_type> (2 * ht / hx, V1, i + 1, j, element_i, G);
                set_off_diag <solver_type> (-2 * ht / hx, V1, i, j, element_i, G);
                set_rhs (GVVn (G, i, j) * (2 + ht / hx * (GVVn (V1, i + 1, j) - GVVn (V1, i, j))) + ht / hx * (GVVn (G, i, j) * GVVn (V1, i, j) - 
                      2.5 * GVVn (G, i + 1, j) * GVVn (V1, i + 1, j) + 2 * GVVn (G, i + 2, j) * GVVn (V1, i + 2, j) - 
                      0.5 * GVVn (G, i + 3, j) * GVVn (V1, i + 3, j) + (2 - GVVn (G, i, j)) * 
                      (GVVn (V1, i, j) - 2.5 * GVVn (V1, i + 1, j) + 2 * GVVn (V1, i + 2, j) - 0.5 * GVVn (V1, i + 3, j))) + 2 * ht * Func_0 (time_step * ht, i * hx, j * hy),
                    G, element_i);
               // set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }

            // right boundary  
            case X_RIGHT:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            // down boundary
            case Y_DOWN:
              {
                // equation for G
                set_diag<solver_type> (2 - ht/hy * GVVn (V2, i, j), G, element_i);
                set_off_diag <solver_type> (ht/hy * GVVn (V2, i, j + 1), G, i, j + 1, element_i, G);
                set_off_diag <solver_type> (2 * ht / hy, V2, i, j + 1, element_i, G);
                set_off_diag <solver_type> (-2 * ht / hy, V2, i, j, element_i, G);
                set_rhs (GVVn (G, i, j) * (2 + ht / hy * (GVVn (V2, i, j + 1) - GVVn (V2, i, j))) + ht / hy * (GVVn (G, i, j) * GVVn (V2, i, j) - 
                      2.5 * GVVn (G, i, j + 1) * GVVn (V2, i, j + 1) + 2 * GVVn (G, i, j + 2) * GVVn (V2, i, j + 2) - 
                      0.5 * GVVn (G, i, j + 3) * GVVn (V2, i, j + 3) + (2 - GVVn (G, i, j)) * 
                      (GVVn (V2, i, j) - 2.5 * GVVn (V2, i, j + 1) + 2 * GVVn (V2, i, j + 2) - 0.5 * GVVn (V2, i, j + 3))) + 2 * ht * Func_0 (time_step * ht, i * hx, j * hy),
                    G, element_i);
                //set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            // up boundary  
            case Y_UP:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            case CORNER_1:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            case CORNER_2:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            case CORNER_3:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            case CORNER_4:
              {
                // equation for G
                set_diag<solver_type> (1, G, element_i);
                set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                // equation for V1
                set_diag<solver_type> (1, V1, element_i);
                set_rhs (0, V1, element_i);

                // equation for V2
                set_diag<solver_type> (1, V2, element_i);
                set_rhs (0, V2, element_i);
                break;
              }
            default:
              assert (false);
              break;
         }
      }
    return 0;
  }

template unsigned int matrix_storage::fill_matrix<solver_own> (unsigned int);
template unsigned int matrix_storage::fill_matrix<solver_eigen> (unsigned int);

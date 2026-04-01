#pragma once

#include "grid.hpp"
#include "parse_command_line.hpp"
#include "func.hpp"
#include <cmath>
#include <vector>
#include <algorithm>


#define G 0
#define V1 1
#define V2 2


class matrix_storage
{
  Grid grid; 

  double* matrix = nullptr;
  unsigned int* I = nullptr;

  double* rhs = nullptr;

  double* GVV = nullptr;
  double* GVVn_ = nullptr;


  unsigned int matrix_size = 0;

  double pp = 0;
  double mu = 0;


  unsigned int get_column_num (unsigned int variable, unsigned int i, unsigned int j)
    {
      unsigned int column = grid.convert_ij_to_element_i (i, j);
      column = column * (variable + 1) + variable;
      return column;
    }

  void set_off_diag (double val, unsigned int variable, unsigned int i /*shifted node */, unsigned int j /* shifted node */ , unsigned int element_i /* current node */)
    {
      unsigned int column = get_column_num (variable, i, j);
      element_i = element_i * (variable + 1) + variable;
      unsigned int l = I[element_i + 1] - I[element_i];
      unsigned int J = I[element_i];
      unsigned int k = 0;
      for (k = 0; k < l; ++k)
        {
          if (column == I[J + k])
            {
              break;
            }
        }
      if (k == l)
        assert (false);
      matrix[J + k] = val;
    }
  void set_diag (double val, unsigned int variable, unsigned int element_i)
    {
      element_i = element_i * (variable + 1) + variable;
      matrix[element_i] = val;
    }
  double GVVn (unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int shifted_element_i = grid.convert_ij_to_element_i (i, j);
        shifted_element_i = shifted_element_i * (variable + 1) + variable;
        return GVVn_[shifted_element_i];
    }

  void set_rhs (double val, unsigned int variable, unsigned int element_i)
    {
        element_i = element_i * (variable + 1) + variable;
        GVVn_[element_i] = val;
    }

  double min_G ()
    {
      unsigned int n_elements = grid.get_n_elements ();
      double res = GVVn (G, 0, 0);
      unsigned int  i = 0;
      unsigned int j = 0;
      for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
        {
          grid.convert_element_i_to_ij (element_i, i, j);
          res = std::min (res, GVVn (G, i, j)); 
        }
      return res;
    }

  double pd (double H, double pp)
    {
      if (std::abs (pp - 1.4) < eps_)
        {
          return pp * std::pow (H, pp - 1);
        }
      return pp;
    }

  unsigned int count_matrix_size (void)
    {
      unsigned int sz = 0;
      matrix_size = sz;
      return sz;
    }

public:

  int init_grid (Parser &parser)
    {
      unsigned int Nx = 0;
      unsigned int Ny = 0;
      if (parser.get ("Nx", Nx) < 0)
        return -1;
      if (parser.get ("Ny", Ny) < 0)
        return -1;
      grid.count_number_of_elements (Nx, Ny);
      grid.check_ij_to_n_elememts_mapping ();

      double hx = 0;
      if (parser.get ("hx", hx) < 0)
        return -1;
      double hy = 0;
      if (parser.get ("hy", hy) < 0)
        return -1;
      double ht = 0;
      if (parser.get ("ht", ht) < 0)
        return -1;
      grid.set_h (hx, hy, ht);

      double pp = 0;
      if (parser.get ("pp", pp) < 0)
        return -1;
      double mu = 0;
      if (parser.get ("mu", mu) < 0)
        return -1;

      return 0;
    }

  int set_diff_params (Parser &parser)
    {
      pp = 0;
      if (parser.get ("pp", pp) < 0)
        return -1;
      mu = 0;
      if (parser.get ("mu", mu) < 0)
        return -1;
      return 0;
    }

  int solve ()
    {
      return 0;
    }


  unsigned int fill_matrix (unsigned int time_step)
    {
      unsigned int n_elements = grid.get_n_elements ();
      double hx = 0;
      double hy = 0;
      double ht = 0;
      grid.get_h (&hx, &hy, &ht);

      for (unsigned int element_i = 0; element_i < n_elements; element_i++)
        {
          unsigned int i = 0;
          unsigned int j = 0;
          grid.convert_element_i_to_ij (element_i, i, j);
          unsigned int border_type = grid.get_bored_type (i, j);

          double min_g = min_G ();
          double mum = mu * std::exp (-min_g);

          switch (border_type)
            {
              // inner domain
              case INNER:
                {
                  // equation for G
                  set_diag (4, G, element_i);
                  set_off_diag (-ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)), G, i - 1, j, element_i); 
                  set_off_diag (-ht / hy * (GVVn (V2, i, j - 1) + GVVn (V2, i, j)), G, i, j - 1, element_i);
                  set_off_diag (ht / hx  * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)), G, i + 1, j, element_i);
                  set_off_diag (ht / hy  * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)), G, i, j + 1, element_i);

                  set_off_diag (-2 * ht / hx, V1, i - 1, j, element_i);
                  set_off_diag (-2 * ht / hy, V2, i, j - 1, element_i);
                  set_off_diag (2 * ht / hx, V1, i + 1, j, element_i);
                  set_off_diag (2 * ht / hy, V2, i, j + 1, element_i);

                  set_rhs (4 * GVVn (G, i, j) + ht * GVVn (G, i, j) * 
                      (1. / hx * (GVVn (V1, i + 1, j) - GVVn (V1, i - 1, j)) + 1. / hy * (GVVn (V2, i, j + 1) - GVVn (V2, i, j - 1))) + 4 * ht * Func_0 (time_step * ht, i * hx, j * hy) , G, element_i);

                  // equation for V1
                  set_diag (6 + 4 * ht * mum * (4 / hx / hx + 3 / hy / hy), V1, element_i);
                  set_off_diag (- (ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)) + 8 * ht * mum / hx / hx), V1, i - 1, j, element_i);
                  set_off_diag (-(3 * ht / 2 / hy * (GVVn (V2, i, j - 1) + GVVn (V2 , i, j)) + 6 * ht * mum / hy / hy), V1, i, j - 1, element_i);
                  set_off_diag (ht / hx * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)) - 8 * ht * mum / hx / hx, V1, i + 1, j, element_i);
                  set_off_diag (3 * ht / 2 / hy * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)) - 6 * ht * mum / hy / hy, V1, i, j + 1, element_i);

                  double H = std::exp (GVVn (G, i, j));
                  set_off_diag (-3 * ht * pd (H, pp) / hx, G, i - 1, j, element_i);
                  set_off_diag (3 * ht * pd (H, pp) / hx, G, i + 1, j, element_i);

                  set_rhs (6 * GVVn (V1, i, j) + 3 * ht / 2 / hy * GVVn (V1, i, j) * (GVVn (V2, i, j + 1) + GVVn (V2, i, j - 1)) + 6 * ht * (mu / H - mum) * (4. / 3 / hx / hx * (GVVn (V1, i + 1, j)
                          - 2 * GVVn (V1, i, j) + GVVn (V1, i - 1, j)) + 1 / hy / hy * (GVVn ( V1, i, j + 1) - 2 * GVVn (V1, i, j) + GVVn (V1, i, j - 1))) + 
                      ht * mu / 2 / H / hx / hy * (GVVn (V2, i + 1, j + 1) - GVVn (V2, i - 1, j + 1) - GVVn (V2, i + 1, j - 1) + GVVn (V2, i - 1, j - 1)) + 
                      6 * ht * Func_1 (time_step * ht, i * hx, j * hy, pp, mu)
                      , V1, element_i);

                  // equation for V2
                  set_diag (6 + 4 * ht * mum * (3 / hx / hx + 4 / hy / hy), V2, element_i);
                  set_off_diag (-(3 * ht / 2 / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)) + 6 * ht * mum / hx / hx), V2, i - 1, j, element_i);
                  set_off_diag (-(ht / hy * (GVVn (V2, i, j - 1) + GVVn (V2, i, j)) + 8 * ht * mum / hy / hy), V2, i, j - 1, element_i);
                  set_off_diag (3 * ht / 2 / hx * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)) - 6 * ht * mum / hx / hx, V2, i + 1, j, element_i);
                  set_off_diag (ht / hy * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)) - 8 * ht * mum / hy / hy, V2, i, j + 1, element_i);

                  set_off_diag (-3 * ht * pd (H, pp) / hy, G, i, j - 1, element_i);
                  set_off_diag (3 * ht * pd (H, pp) / hy, G, i, j + 1, element_i);

                  set_rhs (6 * GVVn (V2, i, j) + 3 * ht / 2 / hx * GVVn (V2, i, j) * (GVVn (V1, i + 1, j) + GVVn (V1, i - 1, j)) + 6 * ht * (mu / H - mum) * (1 / hx / hx *
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
                  set_diag (1, G, element_i);
                  set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                  // equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  // equation for V2
                  set_diag (1, V2, element_i);
                  set_rhs (0, V2, element_i);
                  break;
                }

              // right boundary  
              case X_RIGHT:
                {
                  // equation for G
                  set_diag (1, G, element_i);
                  set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                  // equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  break;
                }
              // down boundary
              case Y_DOWN:
                {
                  // equation for G
                  set_diag (1, G, element_i);
                  set_rhs (g (time_step * ht, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  break;
                }
              // up boundary  
              case Y_UP:
                {
                  //equation for G
                  set_diag (1, G, element_i);
                  set_rhs ( g (time_step, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  break;
                }
              case CORNER_1:
                {
                  //equation for G
                  set_diag (1, G, element_i);
                  set_rhs ( g (time_step, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);            
                }
              case CORNER_2:
                {
                  //equation for G
                  set_diag (1, G, element_i);
                  set_rhs ( g (time_step, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);            
                }
              case CORNER_3:
                {
                  //equation for G
                  set_diag (1, G, element_i);
                  set_rhs ( g (time_step, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);            
                }
              case CORNER_4:
                {
                  //equation for G
                  set_diag (1, G, element_i);
                  set_rhs ( g (time_step, i * hx, j * hy), G, element_i);

                  //equation for V1
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);

                  //equation for V2
                  set_diag (1, V1, element_i);
                  set_rhs (0, V1, element_i);            
                }
              default:
                assert (false);
                break;
           }
        }
      return 0;
    }

  int allocate (void)
    {
      return 0;
    }

  int fill_matrix_pattern (void)
    {
      std::vector<unsigned int> columns;
      columns.reserve (100);

      unsigned int n_elements = grid.get_n_elements ();

      unsigned int filled = 3 * n_elements;
      auto set_colums = [this, &columns, &filled]  (unsigned int column)
      {
        I[column] = filled;
        std::sort (columns.begin (), columns.end ());
        unsigned int prev_col = static_cast <unsigned int> (-1);
        for (auto col : columns)
          {
            // remove dublicates && and diagonal
            if (col != prev_col && col != column)
              {
                I[filled] = col;
                filled++;
                prev_col = col;
              }
          }
        columns.clear ();
      };

      for (unsigned int element_i = 0; element_i < n_elements; element_i++)
        {
          unsigned int i = 0;
          unsigned int j = 0;
          grid.convert_element_i_to_ij (element_i, i, j);
          unsigned int border_type = grid.get_bored_type (i, j);

          unsigned int column = 0;

          switch (border_type)
            {
              case INNER:
                column = get_column_num (G, i, j);
                columns.push_back (column);
                column = get_column_num (G, i + 1, j);
                columns.push_back (column);

                column = get_column_num (G, i, j);
                set_colums (column);

                break;
            }

        }
      assert (filled == matrix_size);
      return 0;
    }
};

#pragma once

#include "grid.hpp"
#include "parse_command_line.hpp"
#include <cmath>


#define G 0
#define V1 1
#define V2 2

static inline double eps = 1e-12;

class matrix_storage
{
  Grid grid; 
  double* matrix = nullptr;
  unsigned int* I = nullptr;
  double* GVV = nullptr;
  double* GVVn_ = nullptr;
  double* rhs = nullptr;

  unsigned int matrix_size = 0;

  double pp = 0;
  double mu = 0;

  void set_matrix_value (double val, unsigned int variable, unsigned int i, unsigned int j)
    {
      unsigned int element_i = grid.convert_ij_to_element_i (i, j);
      element_i = element_i * (variable + 1) + variable;
      unsigned int l = I[element_i + 1] - I[element_i];
      unsigned int J = I[element_i];
      unsigned int k = 0;
      for (k = 0; k < l; ++k)
        {
          if (element_i == I[J + k])
            {
              break;
            }
        }
      if (k == l)
        assert (false);
      matrix[J + k] = val;
    }
  double GVVn (unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int shifted_element_i = grid.convert_ij_to_element_i (i, j);
        shifted_element_i = shifted_element_i * (variable + 1) + variable;
        return GVVn_[shifted_element_i];
    }

  void set_rhs (double val, unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int shifted_element_i = grid.convert_ij_to_element_i (i, j);
        shifted_element_i = shifted_element_i * (variable + 1) + variable;
        GVVn_[shifted_element_i] = val;
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
      if (std::abs (pp - 1.4) < eps)
        {
          return pp * std::pow (H, pp - 1);
        }
      return pp;
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
      grid.set_restrictions (Nx, Ny);
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


  unsigned int fill_matrix ()
    {
      unsigned int n_elements = grid.get_n_elements ();
      for (unsigned int element_i = 0; element_i < n_elements; element_i++)
        {
          unsigned int i = 0;
          unsigned int j = 0;
          grid.convert_element_i_to_ij (element_i, i, j);
          unsigned int border_type = grid.get_bored_type (i, j);

          double hx = 0;
          double hy = 0;
          double ht = 0;
          grid.get_h (&hx, &hy, &ht);

          double mum = std::exp (- min_G ());

          switch (border_type)
            {
              // not border
              case COND_NO:
                {
                  set_matrix_value (4, G, i, j);
                  set_matrix_value (-ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)), G, i - 1, j); 
                  set_matrix_value (-ht / hy * (GVVn (V2, i, j - 1) + GVVn (V2, i, j)), G, i, j - 1);
                  set_matrix_value (ht / hx  * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)), G, i + 1, j);
                  set_matrix_value (ht / hy  * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)), G, i, j + 1);

                  set_matrix_value (-2 * ht / hx, V1, i - 1, j);
                  set_matrix_value (-2 * ht / hy, V2, i, j - 1);
                  set_matrix_value (2 * ht / hx, V1, i + 1, j);
                  set_matrix_value (2 * ht / hy, V2, i, j + 1);

                  set_rhs (4 * GVVn (G, i, j) + ht * GVVn (G, i, j) * 
                      (1. / hx * (GVVn (V1, i + 1, j) - GVVn (V1, i - 1, j)) + 1. / hy * (GVVn (V2, i, j + 1) - GVVn (V2, i, j - 1))), G, i, j);


                  set_matrix_value (6 + 4 * ht * mum * (4 / hx / hx + 3 / hy / hy), V1, i, j);
                  set_matrix_value (- (ht / hx * (GVVn (V1, i - 1, j) + GVVn (V1, i, j)) + 8 * ht * mum / hx / hx), V1, i - 1, j);
                  set_matrix_value (-(3 * ht / 2 / hy * (GVVn (V2, i, j - 1) + GVVn (V2 , i, j)) + 6 * ht * mum / hy / hy), V1, i, j - 1);
                  set_matrix_value (ht / hx * (GVVn (V1, i + 1, j) + GVVn (V1, i, j)) - 8 * ht * mum / hx / hx, V1, i + 1, j);
                  set_matrix_value (3 * ht / 2 / hy * (GVVn (V2, i, j + 1) + GVVn (V2, i, j)) - 6 * ht * mum / hy / hy, V1, i, j + 1);

                  double H = std::exp (GVVn (G, i, j));
                  set_matrix_value (-3 * ht * pd (H, pp) / hx, G, i - 1, j);
                  set_matrix_value (3 * ht * pd (H, pp) / hx, G, i + 1, j);

                  set_rhs (0, V1, i, j);


                  break;
                }
              // zero velocity condition  
              case COND_U0:
                break;
              // velocity condition  
              case COND_U:
                break;
              // normal derivative condition  
              case COND_DU:
                break;
              default:
                assert (false);
                break;
           }
        }
      return 0;
    }

  unsigned int count_matrix_size (void)
    {
      unsigned int n_elements = grid.get_n_elements ();
      unsigned int sz = 0;
      for (unsigned int element_i = 0; element_i < n_elements; element_i++)
        {
          unsigned int i = 0;
          unsigned int j = 0;
          grid.convert_element_i_to_ij (element_i, i, j);
          unsigned int border_type = grid.get_bored_type (i, j);

          switch (border_type)
            {
              // not border
              case COND_NO:
                sz += 1;
                break;
              // zero velocity condition  
              case COND_U0:
                sz += 1;
                break;
              // velocity condition  
              case COND_U:
                sz += 1;
                break;
              // normal derivative condition  
              case COND_DU:
                sz += 1;
                break;
              default:
                assert (false);
                break;
           }
        }
      matrix_size = sz;
      return sz;
    }
};

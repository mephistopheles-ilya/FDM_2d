#pragma once

#include "grid.hpp"
#include "parse_command_line.hpp"
#include "func.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include "func.hpp"
#include "solve.hpp"


#define G 0
#define V1 1
#define V2 2


class matrix_storage
{
  Grid grid; 
  linear_solver solver;

  double* matrix = nullptr;
  unsigned int* I = nullptr;

  double* rhs = nullptr;

  double* GVV_ = nullptr;
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
        unsigned int element_i = grid.convert_ij_to_element_i (i, j);
        element_i = element_i * (variable + 1) + variable;
        return GVVn_[element_i];
    }
  double GVV (unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int element_i = grid.convert_ij_to_element_i (i, j);
        element_i = element_i * (variable + 1) + variable;
        return GVV_[element_i];
    }

  void set_rhs (double val, unsigned int variable, unsigned int element_i)
    {
        element_i = element_i * (variable + 1) + variable;
        rhs[element_i] = val;
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

public:

  void update_prev_solution (void)
  {
    unsigned int element_i = grid.get_n_elements ();
    memcpy (GVVn_, GVV_, 3 * element_i * sizeof (double));
  }

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

  int init_solver ()
  {
    unsigned int n_elements = grid.get_n_elements ();
    int ret;
    ret = solver.init_solver (3 * n_elements);
    return ret;
  }

  int solve ()
    {
      int ret = 0;
      ret = solver.solve (matrix, I, rhs, GVV_);
      return ret;
    }


  unsigned int fill_matrix (unsigned int time_step);

  int allocate (void);

  int fill_matrix_pattern (void);

  ~matrix_storage (void)
    {
      delete[] matrix;
      delete[] I;
      delete[] rhs;
      delete[] GVV_;
      delete[] GVVn_;
    }

  template <unsigned int variable>
  double calculate_C_norm (unsigned int time_step)
  {
    constexpr double (*func) (double, double , double) = variable == G ? g : variable == V1 ? u1 : u2;
    unsigned int n_elements = grid.get_n_elements ();
    unsigned int  i = 0;
    unsigned int j = 0;
    double C_norm = 0;
    double hx = 0, hy = 0, ht = 0;
    grid.get_h (&hx, &hy, &ht);
    for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
      {
        grid.convert_element_i_to_ij (element_i, i, j);
        double calc_val = GVV (variable, i, j);
        double real_val = func (time_step * ht, i * hx, j * hy);
        double diff_abs = fabs (calc_val - real_val);
        C_norm = std::max (C_norm, diff_abs);
      }
    return C_norm;
  }

  template <unsigned int variable>
  double calculate_L2_norm (unsigned int time_step)
  {
    constexpr double (*func) (double, double , double) = variable == G ? g : variable == V1 ? u1 : u2;
    unsigned int n_elements = grid.get_n_elements ();
    unsigned int  i = 0;
    unsigned int j = 0;
    double L2_norm = 0;
    double hx = 0, hy = 0, ht = 0;
    grid.get_h (&hx, &hy, &ht);
    for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
      {
        grid.convert_element_i_to_ij (element_i, i, j);
        double calc_val = GVV (variable, i, j);
        double real_val = func (time_step * ht, i * hx, j * hy);
        double diff = calc_val - real_val;
        diff *= diff;
        L2_norm += diff;
      }
    return sqrt (L2_norm * hx * hy);
  }

  template <unsigned int variable>
  double calculate_W1_norm (unsigned int time_step)
  {
    return time_step;
  }

};

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

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


#define G 0
#define V1 1
#define V2 2
#define VRAS_NUM 3

#define solver_own 0
#define solver_eigen 1


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

  Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_A;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IdentityPreconditioner> eigen_solver;
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::DiagonalPreconditioner<double>> eigen_solver;
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> eigen_solver;
  double *A_values = nullptr;    
  int *A_inner_indices = nullptr;
  int *A_outer_starts = nullptr; 


  unsigned int get_column_num (unsigned int variable, unsigned int i, unsigned int j)
    {
      unsigned int column = grid.convert_ij_to_element_i (i, j);
      column = column * VRAS_NUM + variable;
      return column;
    }

  template <unsigned int solver_type>
  void set_off_diag (double val, unsigned int variable, unsigned int i /*shifted node */, unsigned int j /* shifted node */ , unsigned int element_i /* current node */, unsigned int eq)
    {
      if constexpr (solver_type == solver_own)
        {
          element_i = element_i * VRAS_NUM + eq;
          unsigned int l = I[element_i + 1] - I[element_i];
          unsigned int J = I[element_i];
          unsigned int k = 0;
          unsigned int column = get_column_num (variable, i, j);
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
      if constexpr (solver_type == solver_eigen)
        {
          unsigned int row = element_i * VRAS_NUM + eq;
          unsigned int col = get_column_num (variable, i, j);
          unsigned int idx = static_cast <unsigned int> (A_outer_starts[row]);
          unsigned int idx_end = static_cast <unsigned int> (A_outer_starts[row + 1]);
          for (; idx < idx_end; ++idx)
            {
              if (A_inner_indices[idx] == static_cast<int> (col))
                {
                  break;
                }
            }
          if (idx == idx_end)
            assert (false);
          A_values[idx] = val;
        }
    }
  template <unsigned int solver_type>
  void set_diag (double val, unsigned int variable, unsigned int element_i)
    {
      if constexpr (solver_type == solver_own)
        {
          element_i = element_i * VRAS_NUM + variable;
          matrix[element_i] = val;
        }
      if constexpr (solver_type == solver_eigen)
        {
          unsigned int row = element_i * VRAS_NUM + variable;
          unsigned int col = row;
          unsigned int idx = static_cast <unsigned int> (A_outer_starts[row]);
          unsigned int idx_end = static_cast <unsigned int> (A_outer_starts[row + 1]);
          for (; idx < idx_end; ++idx)
            {
              if (A_inner_indices[idx] == static_cast<int> (col))
                {
                  break;
                }
            }
          if (idx == idx_end)
            assert (false);
          A_values[idx] = val;

        }
    }
  double GVVn (unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int element_i = grid.convert_ij_to_element_i (i, j);
        element_i = element_i * VRAS_NUM + variable;
        return GVVn_[element_i];
    }
  double GVV (unsigned int variable, unsigned int i, unsigned int j)
    {
        unsigned int element_i = grid.convert_ij_to_element_i (i, j);
        element_i = element_i * VRAS_NUM + variable;
        return GVV_[element_i];
    }

  void set_rhs (double val, unsigned int variable, unsigned int element_i)
    {
        element_i = element_i * VRAS_NUM + variable;
        rhs[element_i] = val;
    }

  void set_GVVn (double val, unsigned int variable, unsigned int element_i)
    {
        element_i = element_i * VRAS_NUM + variable;
        GVVn_[element_i] = val;
    }

  void set_GVV (double val, unsigned int variable, unsigned int element_i)
    {
        element_i = element_i * VRAS_NUM + variable;
        GVV_[element_i] = val;
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
      if (std::abs (pp - 1.4) < eq_eps)
        {
          return pp * std::pow (H, pp - 1);
        }
      return pp;
    }

  void init_eigen_from_msr ()
  {
    unsigned int n_elements = grid.get_n_elements ();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve (matrix_size - n_elements);

    for(unsigned int i = 0; i < 3 * n_elements; ++i)
      {
          unsigned int row = i;
          triplets.emplace_back (row, row, 0.0);
          unsigned int l = I[i + 1] - I[i];
          unsigned int J = I[i];
          for(unsigned int j = 0; j < l; ++j)
            {
              unsigned int col = I[J + j];
              triplets.emplace_back (row, col, 0.0);
            }
      }
    eigen_A.resize (3 * n_elements, 3 * n_elements);
    eigen_A.setFromTriplets (triplets.begin(), triplets.end());
    eigen_A.makeCompressed ();

    A_values = eigen_A.valuePtr ();
    A_inner_indices = eigen_A.innerIndexPtr ();
    A_outer_starts = eigen_A.outerIndexPtr ();
  }
public:



  void init_solution (void)
  {
    double hx = 0;
    double hy = 0;
    unsigned int n_elements = grid.get_n_elements ();
    grid.get_h (&hx, &hy, nullptr);
    unsigned int i = 0;
    unsigned int j = 0;
    for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
      {
        grid.convert_element_i_to_ij (element_i, i, j);
        double val = g (0, i * hx, j * hy);
        set_GVVn (val, G, element_i);
        set_GVV (val, G, element_i);

        val = u1 (0, i * hx, j * hy);
        set_GVVn (val, V1, element_i);
        set_GVV (val, V1, element_i);

        val = u2 (0, i * hx, j * hy);
        set_GVVn (val, V2, element_i);
        set_GVV (val, V2, element_i);
      }
  }

  template <unsigned int variable>
  void correct (unsigned int time_step)
  {
    constexpr double (*func) (double, double , double) = variable == G ? g : variable == V1 ? u1 : u2;
    double hx = 0;
    double hy = 0;
    double ht = 0;
    unsigned int n_elements = grid.get_n_elements ();
    grid.get_h (&hx, &hy, &ht);
    unsigned int i = 0;
    unsigned int j = 0;
    for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
      {
        grid.convert_element_i_to_ij (element_i, i, j);
        double val = func (time_step * ht, i * hx, j * hy);
        set_GVV (val, variable, element_i);
      }
  }

  void update_prev_solution (void)
    {
      unsigned int n_elements = grid.get_n_elements ();
      memcpy (GVVn_, GVV_, 3 * n_elements * sizeof (double));
    }

  int init_grid (Parser &parser)
    {
      unsigned int Nx = 0;
      unsigned int Ny = 0;
      if (parser.get ("Nx", Nx) < 0)
        return -1;
      if (parser.get ("Ny", Ny) < 0)
        return -1;
      grid.set_N (Nx, Ny);
      grid.count_number_of_elements ();
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

  int init_solver (Parser &parser)
  {
    unsigned int n_elements = grid.get_n_elements ();
    double eps = 0;
    if (parser.get ("eps", eps) < 0)
      return -1;
    unsigned int maxit = 0; 
    if (parser.get ("maxit", maxit) < 0)
      return -1;
    unsigned int solver_type = 0;
    if (parser.get ("solver", solver_type) < 0)
      return -1;
    if (solver_type == solver_own)
      {
        solver.set_parms (eps, maxit);
        int ret;
        ret = solver.init_solver (3 * n_elements);
        return ret;
      }
    if (solver_type == solver_eigen)
      {
        init_eigen_from_msr ();
        eigen_solver.setTolerance (eps);
        eigen_solver.setMaxIterations (maxit);
      }
    return 0;
  }

  int solve (unsigned int solver_type)
    {
      int ret = 0;
      if (solver_type == solver_own)
        {
          ret = solver.solve (matrix, I, rhs, GVV_);
          return ret;
        }
      else
        {
          eigen_solver.compute (eigen_A);
          if (eigen_solver.info () != Eigen::Success)
            {
              return -1;
            }
          unsigned int n_elements = grid.get_n_elements ();
          Eigen::Map<Eigen::VectorXd> solution_map (GVV_, 3 * n_elements);
          solution_map = eigen_solver.solve (Eigen::Map<Eigen::VectorXd> (rhs, 3 * n_elements));
          if (eigen_solver.info() != Eigen::Success)
            {
              std::cerr << "ERROR: Compute failed!" << std::endl;
              std::cerr << "Reason: " << eigen_solver.info () << std::endl;
              std::cerr << "Last error estimate: " << eigen_solver.error() << std::endl;
              std::cerr << "IT = " << eigen_solver.iterations ();
              return -1;
            }
          return eigen_solver.iterations ();
       }
      return ret;
    }


  template <unsigned int solver_type>
  unsigned int fill_matrix (unsigned int time_step);

  int allocate (unsigned int solver_type);

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
        unsigned int border_type = grid.get_bored_type (i, j);
        if  (border_type == INNER)
          {
            diff *= diff;
          }
        else
          {
            diff *= diff * 0.5;
          }
        L2_norm += diff;
      }
    return sqrt (L2_norm * hx * hy);
  }

  template <unsigned int variable>
  double calculate_W1_norm (unsigned int time_step)
  {
    constexpr double (*func) (double, double , double) = variable == G ? g : variable == V1 ? u1 : u2;
    unsigned int n_elements = grid.get_n_elements ();
    unsigned int  i = 0;
    unsigned int j = 0;
    double W1_norm_x = 0;
    double W1_norm_y = 0;
    double hx = 0, hy = 0, ht = 0;
    grid.get_h (&hx, &hy, &ht);
    for (unsigned int element_i = 0; element_i < n_elements; ++element_i)
      {
        grid.convert_element_i_to_ij (element_i, i, j);
        double calc_val = GVV (variable, i, j);
        double real_val = func (time_step * ht, i * hx, j * hy);
        double diff = calc_val - real_val;
        if (grid.is_active_node (i + 1, j))
          {
            double forward_calc_val = GVV (variable, i + 1, j);
            double forward_real_val = func (time_step * ht, (i + 1) * hx, j * hy);
            double forward_diff = forward_calc_val - forward_real_val;
            double forward_der = (forward_diff - diff) / hx;
            W1_norm_x += forward_der * forward_der;
          }
        if (grid.is_active_node (i, j + 1))
          {
            double forward_calc_val = GVV (variable, i, j + 1);
            double forward_real_val = func (time_step * ht, i * hx, (j + 1) * hy);
            double forward_diff = forward_calc_val - forward_real_val;
            double forward_der = (forward_diff - diff) / hy;
            W1_norm_y += forward_der * forward_der;
          }
      }
    return sqrt ((W1_norm_x + W1_norm_y) * hx * hy);
  }

};

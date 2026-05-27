#pragma once

#include <cassert>
#include <cmath>
#include "func.hpp"


#define X_LEFT   0
#define X_RIGHT  1
#define Y_DOWN   3
#define Y_UP     4
#define INNER    5

class Grid
{
  static inline double eps = 1e-12;

  unsigned int Nx = 0;
  unsigned int Ny = 0;
  unsigned int n_elements = 0;

  unsigned int border1x = 0;
  unsigned int border2x = 0;
  unsigned int border1y = 0;

  unsigned int amount1 = 0;
  unsigned int amount2 = 0;
  unsigned int amount3 = 0;

  double hx = 0;
  double hy = 0;
  double ht = 0;

public:

  unsigned int get_n_elements (void)
    {
      return n_elements;
    }

  void set_h (double hx_, double hy_, double ht_)
    {
      hx = hx_;
      hy = hy_;
      ht = ht_;
    }

  void set_N (unsigned int Nx_, unsigned int Ny_)
    {
      Nx = Nx_;
      Ny = Ny_;
    }

  void get_h (double *hx_, double *hy_, double *ht_)
    {
      if (hx_)
        *hx_ = hx;
      if (hy_)
        *hy_ = hy;
      if (ht_)
        *ht_ = ht;
    }

  void consctruct_ij_to_n_elements_mapping ()
    {
      border1x = std::ceil(1. / hx);
      assert (fabs (border1x * hx - 1.) < 1e-16);
      border2x = std::ceil(2. / hx);
      assert (fabs (border2x * hx - 2.) < 1e-16);
      border1y = std::ceil(1. / hy);
      assert (fabs (border1y * hy - 1.) < 1e-16);

      assert (fabs (Nx * hx - 3) < 1e-16);
      assert (fabs (Ny * hy - 2) < 1e-16);

      for (unsigned int i = border1x; i <= border2x; ++i)
        {
          double u1_on_border1y = u1 (0, hx * i, hy * border1y);
          assert (fabs (u1_on_border1y) < eps);

          double u2_on_border1y = u2 (0, hx * i, hy * border1y);
          assert (fabs (u2_on_border1y) < eps);
        }
      for (unsigned int j = 0; j <= border1y; ++j)
        {
          double u1_on_border1x = u1 (0, hx * border1x, hy * j);
          assert (fabs (u1_on_border1x) < eps);
          double u1_on_border2x = u1 (0, hx * border2x, hy * j);
          assert (fabs (u1_on_border2x) < eps);

          double u2_on_border1x = u2 (0, hx * border1x, hy * j);
          assert (fabs (u2_on_border1x) < eps);
          double u2_on_border2x = u2 (0, hx * border2x, hy * j);
          assert (fabs (u2_on_border2x) < eps);
        }

      amount1 = (border1x + 1) * (Ny + 1);
      amount2 = (border2x - border1x - 1) * (Ny - border1y + 1);
      amount3 = (Nx - border2x + 1) * (Ny + 1);
    }

  bool is_active_node (unsigned int i, unsigned int j)
    {
      if (i > Nx || j > Ny)
        return false;
      if (i > border1x && i < border2x && j < border1y)
        return false;
      return true;
    }
  unsigned int get_bored_type (unsigned int i, unsigned int j)
    {
      assert (is_active_node (i, j));

      if (i == 0)
        return X_LEFT;
      if (i == Nx)
        return X_RIGHT;
      if (j == 0)
        return Y_DOWN;
      if (j == Ny)
        return Y_UP;
      if (i == border1x && j <= border1y)
        return X_RIGHT;
      if (i == border2x && j <= border1y)
        return X_LEFT;
      if (j == border1y && i > border1x && i < border2x)
        return Y_DOWN;

      return INNER;
    }

  unsigned int convert_ij_to_element_i (unsigned int i, unsigned j)
    {
      assert (is_active_node (i, j));
      if (i <= border1x)
        {
          unsigned int element_i = i * (Ny + 1) + j;
          return element_i;
        }
      if (i < border2x)
        {
          unsigned int element_i = amount1 + (i - border1x - 1) * (Ny - border1y + 1) + (j - border1y);
          return element_i;
        }

      unsigned int element_i = amount1 + amount2 + (i - border2x) * (Ny + 1) + j;
      return element_i;

    }
  void convert_element_i_to_ij (unsigned int element_i, unsigned int &i, unsigned int &j)
    {
      if (element_i < amount1)
        {
          i = element_i / (Ny + 1);
          j = element_i % (Ny + 1);
          return;
        }
      if (element_i < amount1 + amount2)
        {
          element_i -= amount1;
          i = element_i / (Ny - border1y + 1);
          i += border1x + 1;
          j = element_i % (Ny - border1y + 1);
          j += border1y;
          return;
        }
      element_i -= amount1 + amount2;
      i = element_i / (Ny + 1);
      i += border2x;
      j = element_i % (Ny + 1);
    }

  unsigned int count_number_of_elements ()
    {
      unsigned int counter = 0;
      for (unsigned int i = 0; i <= Nx; ++i)
        {
          for (unsigned int j = 0; j <= Ny; ++j)
            {
              if (is_active_node (i, j) == false)
                continue;
              counter++;
              if (i == border1x && j == Ny)
                {
                  assert (counter == amount1);
                }
              if (i == (border2x - 1) && j == Ny)
                {
                  assert (counter == amount1 + amount2);
                }
            }
        }
      n_elements = counter;
      assert (counter == amount1 + amount2 + amount3);
      return n_elements;
    }

  void check_ij_to_n_elememts_mapping (void)
    {
      unsigned int counter = 0;
      for (unsigned int i = 0; i <= Nx; ++i)
        {
          for (unsigned int j = 0; j <= Ny; ++j)
            {
              if (is_active_node(i, j) == false)
                  continue;
              unsigned int ij_to_element_i = convert_ij_to_element_i (i, j);
              assert (ij_to_element_i == counter);
              unsigned int i1 = 0, j1 = 0;
              convert_element_i_to_ij (counter, i1, j1);
              assert (i == i1);
              assert (j == j1);
              counter++;
            }
        }
    }


};

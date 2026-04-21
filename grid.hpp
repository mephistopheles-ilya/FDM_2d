#pragma once

#include <cassert>


#define X_LEFT   0
#define X_RIGHT  1
#define Y_DOWN   3
#define Y_UP     4
#define INNER    5

class Grid
{
  unsigned int Nx = 0;
  unsigned int Ny = 0;
  unsigned int n_elements = 0;

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

  bool is_active_node (unsigned int i, unsigned int j)
    {
      if (i > Nx || j > Ny)
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

      return INNER;
    }

  unsigned int convert_ij_to_element_i (unsigned int i, unsigned j)
    {
      assert (is_active_node (i, j));
      unsigned int element_i = i * (Ny + 1) + j;

      return element_i;

    }
  void convert_element_i_to_ij (unsigned int element_i, unsigned int &i, unsigned int &j)
    {
      i = element_i / (Ny + 1);
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
            }
        }
      n_elements = counter;
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

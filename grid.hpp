#pragma once

#include <cassert>


#define X_LEFT   0
#define X_RIGHT  1
#define Y_DOWN   3
#define Y_UP     4
#define INNER    5
#define CORNER_1 6
#define CORNER_2 7
#define CORNER_3 8
#define CORNER_4 9


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
      if (i > Nx || j < Ny)
        return false;
      return true;
    }
  unsigned int get_bored_type (unsigned int i, unsigned int j)
    {
      assert (is_active_node (i, j));

      if (i == 0 && j == 0)
        return CORNER_1;
      if (i == Nx && j == 0)
        return CORNER_2;
      if (i == Nx && j == Ny)
        return CORNER_3;
      if (i == 0 && j == Ny)
        return CORNER_4;

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

      unsigned int full_column = i = (i == 0 ? 0 : i - 1);
      unsigned int element_i = full_column * (Nx + 1) + j;

      return element_i;

    }
  void convert_element_i_to_ij (unsigned int element_i, unsigned int &i, unsigned int &j)
    {
      i = element_i / (Nx + 1);
      j = element_i % (Nx + 1);
    }

  unsigned int count_number_of_elements (unsigned int Nx, unsigned int Ny)
    {
      unsigned int counter = 0;
      for (unsigned int i = 0; i <= Nx; ++i)
        {
          for (unsigned int j = 0; j <= Ny; ++j)
            {
              if (is_active_node (i, j))
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
              if (is_active_node(i, j))
                  continue;
              unsigned int ij_to_n_element = convert_ij_to_element_i (i, j);
              assert (ij_to_n_element == counter);
              unsigned int i1 = 0, j1 = 0;
              convert_element_i_to_ij (counter, i1, j1);
              assert (i == i1);
              assert (j == j1);
            }
        }
    }


};

#pragma once

#include <cassert>


#define COND_NO 1
#define COND_U0 2
#define COND_U 3
#define COND_DU 4



class Grid
{
  unsigned int y1 = 0;
  unsigned int x1 = 0;
  unsigned int x2 = 0;
  unsigned int Nx = 0;
  unsigned int Ny = 0;
  unsigned int n_elements = 0;

public:

  unsigned int get_n_elements (void)
    {
      return n_elements;
    }
  void set_restrictions (unsigned int Nx, unsigned int Ny)
    {
        y1 = Ny / 2;
        x1 = Nx / 3;
        x2 = 2 * Nx / 3;
    }
  bool is_active_node (unsigned int i, unsigned int j)
    {
      if (j < y1 && i < x2 && i > x1)
        return false;
      return true;
    }
  unsigned int get_bored_type (unsigned int i, unsigned int j)
    {
      if (i == 0 && j <= y1)
        return COND_U;
      if (j == 0 && i >= x2)
        return COND_DU;
      if (i == 0 || i == Nx || j == 0 || j == Ny)
        return COND_U0;
      if (j <= y1 && i <= x2 && i >= x1)
        return COND_U0;
      return COND_NO;
    }
  unsigned int convert_ij_to_element_i (unsigned int i, unsigned j)
    {
      if (i <= x1)
        {
          return ((i == 0) ? 0 : (i - 1) * Ny) + j;
        }
      if (i < x2)
        {
          unsigned int n_elements1 = x1 * Ny;
          unsigned int Ny1 = Ny - y1 + 1;
          j  = j - y1;
          i = i - x1;
          return n_elements1 + ((i == 0) ? 0 : (i - 1) * Ny1) + j;
        }
      unsigned int Ny1 = Ny - y1 + 1;
      unsigned int Nx1 = Nx - x1 - (Nx - x2 - 1);
      unsigned int n_elements1 = x1 * Ny + Ny1 * Nx1;
      i = Nx - x2;
      return n_elements1 + ((i == 0) ? 0 : (i - 1) * Ny) + j;
    }
  void convert_element_i_to_ij (unsigned int element_i, unsigned int &i, unsigned int &j)
    {
      unsigned int n_elements1 = x1 * Ny;
      if (element_i <= n_elements1)
        {
          j = element_i % (Ny + 1);
          i = element_i / Nx - (j == Ny ? 1 : 0);
        }
      unsigned int Ny1 = Ny - y1 + 1;
      unsigned int Nx1 = Nx - x1 - (Nx - x2 - 1);
      n_elements1 = x1 * Ny + Ny1 * Nx1;
      if (element_i <= n_elements1)
        {
          element_i -= x1 * Ny;
          j = y1 + element_i % (Ny1 + 1);
          i = x1 + element_i / Ny1 - (j == Ny ? 1 : 0);
        }
      element_i -= n_elements1;
      j = element_i % (Ny + 1);
      i = x2 + element_i / Ny - (j == Ny ? 1 : 0);
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
              unsigned int ij_to_n_element = convert_ij_to_element_i(i, j);
              assert (ij_to_n_element == counter);
              unsigned int i1 = 0, j1 = 0;
              convert_element_i_to_ij(counter, i1, j1);
              assert (i == i1);
              assert (j == j1);
            }
        }
    }


};

#pragma once


#define COND_0 1
#define COND_1 2
#define COND_2 3
#define COND_3 4



class Grid
{
  unsigned int n_elements = 0;
  unsigned int restrict_y = 0;
  unsigned int restrict_x1 = 0;
  unsigned int restrict_x2 = 0;
  unsigned int n_elements_1 = 0;
  unsigned int n_elements_2 = 0;
  unsigned int Nx = 0;
  unsigned int Nx_mid = 0;
  unsigned int Ny = 0;
  unsigned int Ny_up = 0;

  unsigned int empty_element = static_cast<unsigned int > (-1);

public:

  unsigned int get_n_elements (void)
    {
      return n_elements;
    }
  void set_restrictions (unsigned int Nx, unsigned int Ny)
    {
        restrict_y = Ny / 2;
        Ny_up = Ny - restrict_y + 1;
        restrict_x1 = Nx / 3;
        restrict_x2 = 2 * Nx / 3;
        Nx_mid = Nx - restrict_x1 + Nx - restrict_x2;
        n_elements_1 = restrict_x1 * Ny;
        n_elements_2 = restrict_x1 * Ny + Nx_mid * Ny_up; 
    }
  bool is_active_node (unsigned int i, unsigned int j)
    {
      if (j < restrict_y && i < restrict_x2 && i > restrict_x1)
        return false;
      return true;
    }
  unsigned int get_bored_type (unsigned int i, unsigned int j)
    {
      if (i == 0 && j <= restrict_y)
        return COND_2;
      if (j == 0 && i >= restrict_x2)
        return COND_3;
      if (i == 0 || i == Nx || j == 0 || j == Ny)
        return COND_1;
      if (j <= restrict_y && i <= restrict_x2 && i >= restrict_x1)
        return COND_1;
      return COND_0;
    }
  unsigned int convert_ij_to_element_i (unsigned int i, unsigned j)
    {
      if (!is_active_node (i, j))
        return empty_element;
      return 0;
    }
  void convert_element_i_to_ij (unsigned int element_i, unsigned int &i, unsigned int &j)
    {
      if (element_i <= n_elements_1)
        {
          j = element_i % (Ny + 1);
          i = (element_i + 1) / Nx;
        }
      if (element_i < n_elements_2)
        {
          element_i -= n_elements_1;
          j = restrict_y + element_i % (Ny_up + 1);
          i = restrict_x1 + (element_i + 1);
        }
      element_i -= n_elements_2;
      j = element_i % (Ny + 1);
      i = restrict_x2 + (element_i + 1);
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
};

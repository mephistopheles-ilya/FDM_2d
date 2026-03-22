#pragma once

#include "grid.hpp"
#include "parse_command_line.hpp"

class matrix_storage
{
  Grid grid; 
  double* matrix = nullptr;
  double* matrix_pattern = nullptr;
  double* GVV = nullptr;
  double* rhs = nullptr;

  unsigned int matrix_size = 0;

public:

  void init_grid (Parser &parser)
    {
      unsigned int Nx = 0;
      unsigned int Ny = 0;
      parser.get("Nx", Nx);
      parser.get("Ny", Ny);
      grid.set_restrictions (Nx, Ny);
      grid.count_number_of_elements (Nx, Ny);
      grid.check_ij_to_n_elememts_mapping ();
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

          switch (border_type)
            {
              // not border
              case COND_NO:
                break;
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

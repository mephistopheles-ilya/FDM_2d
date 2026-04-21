#include "matrix_storage.hpp"

int matrix_storage::allocate (unsigned int solver_type)
  {
    unsigned int n_elements = grid.get_n_elements ();
    unsigned int sz = 3 * n_elements + 1;
    for (unsigned int element_i = 0; element_i < n_elements; element_i++)
      {
        unsigned int i = 0;
        unsigned int j = 0;
        grid.convert_element_i_to_ij (element_i, i, j);
        unsigned int border_type = grid.get_bored_type (i, j);
        switch (border_type)
          {
            // inner domain
            case INNER:
              {
                // equation for G
                sz += 8;

                // equation for V1
                sz += 6;

                // equation for V2
                sz += 6;

                break;
              }

            // left boundary
            case X_LEFT:
              {
                sz += 3;
                break;
              }

            // right boundary  
            case X_RIGHT:
              {
                sz += 3;
                break;
              }
            // down boundary
            case Y_DOWN:
              {
                sz += 3;
                break;
              }
            // up boundary  
            case Y_UP:
              {
                break;
              }
            case CORNER_1:
              {
                break;
              }
            case CORNER_2:
              {
                break;
              }
            case CORNER_3:
              {
                break;
              }
            case CORNER_4:
              {
                break;
              }
            default:
              assert (false);
              break;
         }
      }
    matrix_size = sz;

    if (solver_type == solver_own)
      matrix = new double [matrix_size];
    I = new unsigned int [matrix_size];
    rhs = new double [3 * n_elements];

    GVV_ = new double [3 * n_elements];
    GVVn_ = new double [3 * n_elements];
    if (I == nullptr || rhs == nullptr || GVV_ == nullptr || GVVn_ == nullptr)
      {
        return -1;
      }
    return 0;
  }

int matrix_storage::fill_matrix_pattern (void)
  {
    std::vector<unsigned int> one_row;
    one_row.reserve (100);

    unsigned int n_elements = grid.get_n_elements ();

    unsigned int filled = 3 * n_elements + 1;
    auto finilize_one_row = [this, &one_row, &filled]  (unsigned int variable, unsigned int i, unsigned int j)
    {
      unsigned int column = get_column_num (variable, i, j);
      I[column] = filled;
      std::sort (one_row.begin (), one_row.end ());
      unsigned int prev_col = static_cast <unsigned int> (-1);
      for (auto col : one_row)
        {
          // must be no diagonal
          assert (col != column);
          // must be no dublicates
          assert (col != prev_col);
          I[filled] = col;
          filled++;
          prev_col = col;
        }
      one_row.clear ();
    };

    auto add_column = [this, &one_row] (unsigned int variable, unsigned i, unsigned j)
    {
      unsigned int column = get_column_num (variable, i, j);
      one_row.push_back (column);
    };

    for (unsigned int element_i = 0; element_i < n_elements; element_i++)
      {
        unsigned int i = 0;
        unsigned int j = 0;
        grid.convert_element_i_to_ij (element_i, i, j);
        unsigned int border_type = grid.get_bored_type (i, j);

        switch (border_type)
          {
            // inner domain
            case INNER:
              {
                // equation for G
                add_column (G, i - 1, j);
                add_column (G, i, j - 1);
                add_column (G, i + 1, j);
                add_column (G, i, j + 1);
                add_column (V1, i - 1, j);
                add_column (V2, i, j - 1);
                add_column (V1, i + 1, j);
                add_column (V2, i, j + 1);

                finilize_one_row (G, i, j);

                // equation for V1
                add_column (V1, i - 1, j);
                add_column (V1, i, j - 1);
                add_column (V1, i + 1, j);
                add_column (V1, i, j + 1);
                add_column (G, i - 1, j);
                add_column (G, i + 1, j);

                finilize_one_row (V1, i, j);

                // equation for V2
                add_column (V2, i - 1, j);
                add_column (V2, i, j - 1);
                add_column (V2, i + 1, j);
                add_column (V2, i, j + 1);
                add_column (G, i, j - 1);
                add_column (G, i, j + 1);

                finilize_one_row(V2, i, j);

                break;
              }

            // left boundary
            case X_LEFT:
              {
                add_column (G, i + 1, j);
                add_column (V1, i + 1, j);
                add_column (V1, i, j);
                finilize_one_row (G, i, j);

                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }

            // right boundary  
            case X_RIGHT:
              {
                add_column (G, i - 1, j);
                add_column (V1, i, j);
                add_column (V1, i - 1, j);
                finilize_one_row (G, i, j);

                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            // down boundary
            case Y_DOWN:
              {
                add_column (G, i, j + 1);
                add_column (V2, i, j + 1);
                add_column (V2, i, j);
                finilize_one_row (G, i, j);

                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            // up boundary  
            case Y_UP:
              {
                finilize_one_row (G, i, j);
                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            case CORNER_1:
              {
                finilize_one_row (G, i, j);
                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            case CORNER_2:
              {
                finilize_one_row (G, i, j);
                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            case CORNER_3:
              {
                finilize_one_row (G, i, j);
                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            case CORNER_4:
              {
                finilize_one_row (G, i, j);
                finilize_one_row (V1, i, j);
                finilize_one_row (V2, i, j);
                break;
              }
            default:
              assert (false);
              break;
         }
      }
    assert (filled == matrix_size);
    I[3 * n_elements] = filled;
    return 0;
  }


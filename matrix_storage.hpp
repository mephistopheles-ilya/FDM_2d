#pragma once

#include "grid.hpp"

class matrix_storage
{
  Grid grid; 
  double* matrix = nullptr;
  double* matrix_pattern = nullptr;
  double* GVV = nullptr;
  double* rhs = nullptr;
};

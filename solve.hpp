#pragma once

#include <cmath>
#include <stdio.h>

inline double min_division = 1e-32;


class linear_solver
{
  unsigned int n = 0;
  double eps = 0;
  unsigned int maxit = 0;
  
  double *z = nullptr;
  double *r = nullptr;
  double *p = nullptr;
  double *u = nullptr;
  double *q = nullptr;
  double *Avec = nullptr;


  void mat_mult_vec (const double *A, const unsigned int *I,  double *vec, double *res)
  {
    for(unsigned int i = 0; i < n; ++i)
      {
          double s = A[i] * vec[i];
          unsigned int l = I[i + 1] - I[i];
          unsigned int J = I[i];
          for(unsigned int j = 0; j < l; ++j)
            {
                s += A[J + j] * vec[I[J + j]];
            }
          res[i] = s;
      }
  }

  double dot (const double *vec1, const double *vec2)
  {
    double s = 0;
    for (unsigned int i = 0; i < n; ++i)
      {
        s += vec1[i] * vec2[i];
      }
    return s;
  }

  void mult_sub_vec (double *res, const double *vec1, double w, const double *vec2)
  {
    for (unsigned int i = 0; i < n; ++i)
      {
        res[i] = vec1[i] - w * vec2[i];
      }
  }

  void copy_vec (double *dest, double *source)
  {
    for (unsigned int i = 0; i < n; ++i)
      {
        dest[i] = source[i];
      }
  }

  void apply_precond_J (const double *A, double *vec)
  {
    for (unsigned int i = 0; i < n; ++i)
      {
        vec[i] /= A[i];
      }
  }


public:

  int init_solver (unsigned int n_)
  {
    n = n_;
    z = new double[n];
    r = new double[n];
    p = new double[n];
    u = new double[n];
    q = new double[n];
    Avec = new double[n];
    if (z == nullptr || r == nullptr || p == nullptr || u == nullptr || q == nullptr || Avec == nullptr)
      return -1;
    return 0;
  }

  void set_parms (double eps_, unsigned int maxit_)
  {
    eps = eps_;
    maxit = maxit_;
  }

  int solve (const double *A, const unsigned int *I, const double *b, double *x)
  {
    mat_mult_vec (A, I, x, Avec);
    mult_sub_vec (r, b, 1, Avec);
    apply_precond_J (A, r);
    copy_vec (p, r);
    copy_vec (u, r);
    copy_vec (z, r);

    double rhs_norm = sqrt (dot (b, b));

    unsigned int it = 0;
    for (it = 0; it < maxit; ++it)
      {
        double rz = dot (r, z);
        if (fabs (rz) < min_division)
          return -1;

        mat_mult_vec (A, I, p, Avec);
        apply_precond_J (A, Avec);
        double Apz = dot (Avec, z);
        if (fabs (Apz) < min_division)
          return -1;

        double alpha = rz / Apz;
        mult_sub_vec (q, u, alpha, Avec);
        mult_sub_vec (u, u, -1, q);
        mult_sub_vec (x, x, -alpha, u);
        mat_mult_vec (A, I, u, Avec);
        apply_precond_J (A, Avec);
        mult_sub_vec (r, r, alpha, Avec);
        double r_norm = sqrt (dot (r, r));
        if (r_norm < eps * rhs_norm)
          return it;
        double rrz = dot (r, z);
        double betta = rrz / rz;
        mult_sub_vec (u, r, -betta, q);
        mult_sub_vec (Avec, q, -betta, p);
        mult_sub_vec (p, u, -betta, Avec);
      }
    return maxit;
  }

  ~linear_solver (void)
  {
    delete[] z;
    delete[] r;
    delete[] p;
    delete[] u;
    delete[] q;
    delete[] Avec;
  }

};

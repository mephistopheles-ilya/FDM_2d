#pragma once

#include <cmath>
#include <stdio.h>
#include <immintrin.h>

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


  void mat_mult_vec (const double *A, const unsigned int *I,  const double *vec, double *res)
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
#if 0
  double dot (const double *vec1, const double *vec2)
  {
    double sum = 0;
    double res256[] = {0, 0, 0, 0};
    unsigned int i = 0;
    __m256d sum0 = _mm256_setzero_pd ();
    __m256d sum1 = _mm256_setzero_pd ();
    __m256d sum2 = _mm256_setzero_pd ();
    __m256d sum3 = _mm256_setzero_pd ();
    for (i = 0; i <= n - 16; i += 16)
      {
        __m256d v1_0 = _mm256_loadu_pd (vec1 + i);
        __m256d v2_0 = _mm256_loadu_pd (vec2 + i);
        __m256d prod0 = _mm256_mul_pd (v1_0, v2_0);
        sum0 = _mm256_add_pd (sum0, prod0);

        __m256d v1_1 = _mm256_loadu_pd (vec1 + i + 4);
        __m256d v2_1 = _mm256_loadu_pd (vec2 + i + 4);
        __m256d prod1 = _mm256_mul_pd (v1_1, v2_1);
        sum1 = _mm256_add_pd (sum1, prod1);

        __m256d v1_2 = _mm256_loadu_pd (vec1 + i + 8);
        __m256d v2_2 = _mm256_loadu_pd (vec2 + i + 8);
        __m256d prod2 = _mm256_mul_pd (v1_2, v2_2);
        sum2 = _mm256_add_pd (sum2, prod2);

        __m256d v1_3 = _mm256_loadu_pd (vec1 + i + 12);
        __m256d v2_3 = _mm256_loadu_pd (vec2 + i + 12);
        __m256d prod3 = _mm256_mul_pd (v1_3, v2_3);
        sum3 = _mm256_add_pd (sum3, prod3);
      }

    __m256d total1 = _mm256_add_pd (sum0, sum1);
    __m256d total2 = _mm256_add_pd (sum2, sum3);
    __m256d total = _mm256_add_pd (total1, total2);

    for (;i < n; ++i)
      {
        sum += vec1[i] * vec2[i];
      }
    _mm256_storeu_pd (res256, total);

    return sum + res256[0] + res256[1] + res256[2] + res256[3];
}
#endif
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
  void linear_combination (double *u, double *p, const double *r, const double *q, double betta)
  {
    double b2 = betta * betta;
    for (unsigned int i = 0; i < n; ++i)
      {
        double bq = betta * q[i];
        double res = r[i] + bq;
        u[i] = res;
        p[i] = res + bq + b2 * p[i];
      }
  }
 void add_vec (double *res, const double *vec)
 {
   for (unsigned int i = 0; i < n; ++i)
    {
      res[i] += vec[i];
    }
 }
 void add_vec (double *res, const double *vec, double alpha)
 {
   for (unsigned int i = 0; i < n; ++i)
    {
      res[i] += alpha * vec[i];
    }
 }

  void copy_vec (double *dest, const double *source)
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
    copy_vec (p, r);
    copy_vec (u, r);
    copy_vec (z, r);

    double rhs_norm = sqrt (dot (b, b));

    double r_norm = 0;
    unsigned int it = 0;
    double rz = dot (r, z);
    double rrz = 0;
    for (it = 0; it < maxit; ++it)
      {
        if (fabs (rz) < min_division)
          return -1;

        mat_mult_vec (A, I, p, Avec);
        double Apz = dot (Avec, z);
        if (fabs (Apz) < min_division)
          return -1;

        double alpha = rz / Apz;
        mult_sub_vec (q, u, alpha, Avec);
        add_vec (u, q);
        add_vec (x, u, alpha);
        mat_mult_vec (A, I, u, Avec);
        mult_sub_vec (r, r, alpha, Avec);
        r_norm = sqrt (dot (r, r));
        if (r_norm < eps * rhs_norm)
          return it;
        rrz = dot (r, z);
        double betta = rrz / rz;
        linear_combination (u, p, r, q, betta);
        rz = rrz;
      }
    printf("last eps = %e, it = %u\n", r_norm / rhs_norm, it);
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

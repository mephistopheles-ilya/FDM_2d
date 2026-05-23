#pragma once

#include <cmath>

static inline double eq_eps = 1e-12;

inline double u1 (double t, double x, double y)
{
  return sin (3 * M_PI * x) * sin (2 * M_PI * y) * exp (t);
}

inline double u2 (double t, double x, double y)
{
  return sin (3 * M_PI * x) * sin (2 * M_PI * y) * exp (-t);
}

inline double rho (double t, double x, double y)
{
  return (cos (2 * M_PI * x) + 1.5) * (sin (2 * M_PI * y) + 1.5) * exp (t);
}

inline double g (double t, double x, double y)
{
  return log (cos (2 * M_PI * x) + 1.5) + log (sin (2 * M_PI * y) + 1.5) + t;
}

inline double dg_dt (double /*t*/, double /*x*/, double /*y*/)
{
  return 1;
}

inline double dg_dx (double /*t*/, double x, double /*y*/)
{
  double num = -2 * M_PI * sin (2 * M_PI * x);
  double den = cos (2 * M_PI * x) + 1.5;
  return num / den;
}

inline double dg_dy (double /*t*/, double /*x*/, double y)
{
  double num = 2 * M_PI * cos (2 * M_PI * y);
  double den = sin (2 * M_PI * y) + 1.5;
  return num / den;
}

inline double du1_dx (double t, double x, double y)
{
  return 3 * M_PI * cos (3 * M_PI * x) * sin (2 * M_PI * y) * exp (t);
}

inline double du1_dy (double t, double x, double y)
{
  return sin (3 * M_PI * x) * 2 * M_PI * cos (2 * M_PI * y) * exp (t);
}

inline double du2_dx (double t, double x, double y)
{
  return 3 * M_PI * cos (3 * M_PI * x) * sin (2 * M_PI * y) * exp (-t);
}

inline double du2_dy (double t, double x, double y)
{
  return sin (3 * M_PI * x) * 2 * M_PI * cos (2 * M_PI * y) * exp (-t);
}

inline double du1_dt (double t, double x, double y)
{
  return u1 (t, x, y);
}

inline double du2_dt (double t, double x, double y)
{
  return - u2 (t, x, y);
}


inline double du1g_dx (double t, double x, double y)
{
  return du1_dx (t, x, y) * g (t, x, y) + u1 (t, x, y) * dg_dx (t, x, y);
}

inline double du2g_dy (double t, double x, double y)
{
  return du2_dy (t, x, y) * g (t, x, y) + u2 (t, x, y) * dg_dy (t, x, y);
}

inline double du1u1_dx (double t, double x, double y)
{
  return 2 * u1 (t, x, y) * du1_dx (t, x, y);
}

inline double du1u2_dy (double t, double x, double y)
{
  return u1 (t, x, y) * du2_dy (t, x, y) + du1_dy (t, x, y) * u2 (t, x, y);
}

inline double dp_drho (double t, double x, double y, double pp)
{
  if (fabs (pp - 1.4) < eq_eps)
    return pp * pow (rho (t, x, y), pp - 1);
  return pp;
}

inline double ddu1_dxdx (double t, double x, double y)
{
  return -9 * M_PI * M_PI * u1 (t, x, y);
}

inline double ddu1_dydy (double t, double x, double y)
{
  return -4 * M_PI * M_PI * u1 (t, x, y);
}

inline double ddu2_dxdy (double t, double x, double y)
{
  return 6 * M_PI * M_PI * cos (3 * M_PI * x) * cos (2 * M_PI * y) * exp (-t);
}

inline double du2u2_dy (double t, double x, double y)
{
  return 2 * u2 (t, x, y) * du2_dy (t, x, y);
}

inline double du1u2_dx (double t, double x, double y)
{
  return u1 (t, x, y) * du2_dx (t, x, y) + du1_dx (t, x, y) * u2 (t, x, y); 
}

inline double ddu2_dydy (double t, double x, double y)
{
  return -4 * M_PI * M_PI * u2 (t, x, y);
}

inline double ddu2_dxdx (double t, double x, double y)
{
  return -9 * M_PI * M_PI * u2 (t, x, y);
}

inline double ddu1_dxdy (double t, double x, double y)
{
  return 6 * M_PI * M_PI * cos (3 * M_PI * x) * cos (2 * M_PI * y) * exp (t);
}

inline double Func_0 (double t, double x, double y)
  {
    return dg_dt (t, x, y) + 
      0.5 * (u1 (t, x, y) * dg_dx (t, x, y) + du1g_dx (t, x, y) + (2 - g (t, x, y)) * du1_dx (t, x, y)) + 
      0.5 * (u2 (t, x, y) * dg_dy (t, x, y) + du2g_dy (t, x, y) + (2 - g (t, x, y)) * du2_dy (t, x, y)); 
  }
inline double Func_1 (double t, double x, double y, double pp, double mu)
  {
    return du1_dt (t, x, y) + 
      1./3. * (u1 (t, x, y) * du1_dx (t, x, y) + du1u1_dx (t, x, y)) + 
      1./2. * (u2 (t, x, y) * du1_dy (t, x, y) + du1u2_dy (t, x, y) - u1 (t, x, y) * du2_dy (t, x, y)) + 
      dp_drho (t, x, y, pp) * dg_dx (t, x, y) - 
      mu / rho (t, x, y) * (4./3. * ddu1_dxdx (t, x, y) + ddu1_dydy (t, x, y) + 1./3. * ddu2_dxdy (t, x, y));
  }
inline double Func_2 (double t, double x, double y, double pp, double mu)
  {
    return du2_dt (t, x, y) + 
      1./3. * (u2 (t, x, y) * du2_dy (t, x, y) + du2u2_dy (t, x, y)) + 
      1./2. * (u1 (t, x, y) * du2_dx (t, x, y) + du1u2_dx (t, x, y) - u2 (t, x, y) * du1_dx (t, x, y)) + 
      dp_drho (t, x, y, pp) * dg_dy (t, x, y) -
      mu / rho (t, x, y) * (4./3. * ddu2_dydy (t, x, y) + ddu2_dxdx (t, x, y) + 1./3. * ddu1_dxdy (t, x, y));
  }


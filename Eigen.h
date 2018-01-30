#ifndef _H_EIGEN
#define _H_EIGEN

  long eigen (double[3][3], double[3][3], double[3],int checkOrth=0);
  double trace (double[3][3]);
  double trace2 (double[3][3]);
  double det (double[3][3]);
  long FindCubicRoots(const double[4], double[3]);
  long FindQuadraticRoots(const double, const double, double[3]);
  long NullSpace(const double *, double *, double, long);
  int checkUnitaryOthoganal(double e[3][3],int &factor);
#ifdef SIM
  double dotProd(double [],double []);
  void crossProd(double [],double [],double []);
  void diffVt(double [],double [],double []);
  int normVt(double [],double[]);
#endif
#endif


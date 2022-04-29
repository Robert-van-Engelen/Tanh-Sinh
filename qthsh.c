// qthsh: Tanh-Sinh quadrature formula
// https://www.genivia.com/qthsh.html
// Dr. Robert A. van Engelen

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// compile with -DFAST to accellerate with modest descrease in accuracy
#ifdef FAST
#define FUDGE1 160
#define FUDGE2 16
#else
#define FUDGE1 10
#define FUDGE2 1
#endif

// integrate function f, range a..b, max levels n (2 to 7, 6 is recommended), relative error tolerance eps, estimated relative error err
double qthsh(double (*f)(double), double a, double b, int n, double eps, double *err) {
  const double tol = FUDGE1*eps;
  double c = (a+b)/2;
  double d = (b-a)/2;
  double s = f(c);
  double p, e, v, h = 2;
  int k = 0;
  do {
    double p = 0, q, fp = 0, fm = 0, t, eh;
    h /= 2;
    eh = exp(h);
    t = eh;
    if (k > 0)
      eh *= eh;
    do {
      double u = exp(1/t-t);      // = exp(-2*sinh(j*h)) = 1/exp(sinh(j*h))^2
      double r = 2*u/(1+u);       // = 1 - tanh(sinh(j*h))
      double w = (t+1/t)*r/(1+u); // = cosh(j*h)/cosh(sinh(j*h))^2
      double x = d*r;
      if (a+x > a) {              // if too close to a then reuse previous fp
        double y = f(a+x);
        if (isfinite(y))
          fp = y;                 // if f(x) is finite, add to local sum
      }
      if (b-x < b) {              // if too close to b then reuse previous fm
        double y = f(b-x);
        if (isfinite(y))
          fm = y;                 // if f(x) is finite, add to local sum
      }
      q = w*(fp+fm);
      p += q;
      t *= eh;
    } while (fabs(q) > eps*fabs(p));
    v = s-p;
    s += p;
    ++k;
  } while (fabs(v) > tol*fabs(s) && k <= n);
  // if the estimated relative error is desired, then return it
  if (err != NULL)
    *err = fabs(v)/(FUDGE2*fabs(s)+eps);
  // result with estimated relative error err
  return d*s*h;
}

// example function to integrate, exact integral x=0..1 is 1
double f(double x) { return acos(x); }

int main() {
  printf("integrate(acos(x), x=0..1) = %.15g\n", qthsh(f, 0, 1, 6, 1e-9, NULL));
}

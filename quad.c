// quad: Tanh-Sinh, Sinh-Sinh and Exp-Sinh quadrature formula
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

#define sign(x) (((x)>0)-((x)<0))

double exp_sinh_opt_d(double (*f)(double), double a, double eps, double d) {
  int ev = 2;
  // const double base = 2; // 2 or 3 or exp(1) for example
  double h2 = f(a + d/2) - f(a + d*2)*4;
  int i = 1, j = 32;     // j=32 is optimal to search for r
  if (isfinite(h2) && fabs(h2) > 1e-5) { // if |h2| > 2^-16
    double r, fl, fr, h, s = 0, lfl, lfr, lr = 2;
    do {                 // find max j such that fl and fr are finite
      j /= 2;
      r = 1 << (i + j);
      fl = f(a + d/r);
      fr = f(a + d*r)*r*r;
      ev += 2;
      h = fl - fr;
    } while (j > 1 && !isfinite(h));
    if (j > 1 && isfinite(h) && sign(h) != sign(h2)) {
      lfl = fl;          // last fl=f(a+d/r)
      lfr = fr;          // last fr=f(a+d*r)*r*r
      do {               // bisect in 4 iterations
        j /= 2;
        r = 1 << (i + j);
        fl = f(a + d/r);
        fr = f(a + d*r)*r*r;
        ev += 2;
        h = fl - fr;
        if (isfinite(h)) {
          s += fabs(h);  // sum |h| to remove noisy cases
          if (sign(h) == sign(h2)) {
            i += j;      // search right half
          }
          else {         // search left half
            lfl = fl;    // record last fl=f(a+d/r)
            lfr = fr;    // record last fr=f(a+d*r)*r*r
            lr = r;      // record last r
          }
        }
      } while (j > 1);
      if (s > eps) {     // if sum of |h| > eps
        h = lfl - lfr;   // use last fl and fr before the sign change
        r = lr;          // use last r before the sign change
        if (h != 0)      // if last difference was nonzero, back up r by one step
          r /= 2;
        if (fabs(lfl) < fabs(lfr))
          d /= r;        // move d closer to the finite endpoint
        else
          d *= r;        // move d closer to the infinite endpoint
      }
    }
  }
  return d;
}

// integrate function f, range a..b, max levels n (2 to 7, 6 is recommended), relative error tolerance eps, estimated relative error err
double quad(double (*f)(double), double a, double b, int n, double eps, double *err) {
  const double tol = FUDGE1*eps;
  double c = 0, d = 1, s, sign = 1, e, v, h = 2;
  int k = 0, mode = 0; // Tanh-Sinh = 0, Exp-Sinh = 1, Sinh-Sinh = 2
  if (b < a) { // swap bounds
    v = b;
    b = a;
    a = v;
    sign = -1;
  }
  if (isfinite(a) && isfinite(b)) {
    c = (a+b)/2;
    d = (b-a)/2;
    v = c;
  }
  else if (isfinite(a)) {
    mode = 1; // Exp-Sinh
    d = exp_sinh_opt_d(f, a, eps, d);
    c = a;
    v = a+d;
  }
  else if (isfinite(b)) {
    mode = 1; // Exp-Sinh
    d = exp_sinh_opt_d(f, b, eps, -d);
    sign = -sign;
    c = b;
    v = b+d;
  }
  else {
    mode = 2; // Sinh-Sinh
    v = 0;
  }
  s = f(v);
  do {
    double p = 0, q, fp = 0, fm = 0, t, eh;
    h /= 2;
    t = eh = exp(h);
    if (k > 0)
      eh *= eh;
    if (mode == 0) {                // Tanh-Sinh
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
    }
    else {
      t /= 2;
      do {
        double r = exp(t-.25/t); // = exp(sinh(j*h))
        double x, y, w = r;
        q = 0;
        if (mode == 1) {         // Exp-Sinh
          x = c + d/r;
          if (x == c)            // if x hit the finite endpoint then break
            break;
          y = f(x);
          if (isfinite(y))       // if f(x) is finite, add to local sum
            q += y/w;
        }
        else {                   // Sinh-Sinh
          r = (r-1/r)/2;         // = sinh(sinh(j*h))
          w = (w+1/w)/2;         // = cosh(sinh(j*h))
          x = c - d*r;
          y = f(x);
          if (isfinite(y))       // if f(x) is finite, add to local sum
            q += y*w;
        }
        x = c + d*r;
        y = f(x);
        if (isfinite(y))         // if f(x) is finite, add to local sum
          q += y*w;
        q *= t+.25/t;            // q *= cosh(j*h)
        p += q;
        t *= eh;
      } while (fabs(q) > eps*fabs(p));
    }
    v = s-p;
    s += p;
    ++k;
  } while (fabs(v) > tol*fabs(s) && k <= n);
  // if the estimated relative error is desired, then return it
  if (err != NULL)
    *err = fabs(v)/(FUDGE2*fabs(s)+eps);
  // result with estimated relative error err
  return sign*d*s*h;
}

// example functions to integrate, exact integrals are 1, 5 and 2
double f1(double x) { return acos(x); }
double f2(double x) { return exp(-x/5); }
double f3(double x) { return pow(cosh(x),-2); }

int main() {
  printf("integrate(acos(x), x=0..1) = %.15g\n", quad(f1, 0, 1, 6, 1e-9, NULL));
  printf("integrate(exp(-x/5), x=0..+inf) = %.15g\n", quad(f2, 0, INFINITY, 6, 1e-9, NULL));
  printf("integrate(1/cosh(x)^2, x=-inf..+inf) = %.15g\n", quad(f3, -INFINITY, INFINITY, 6, 1e-9, NULL));
}

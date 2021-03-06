# Fast double exponential quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh formulas

This work significantly improves the *Tanh-Sinh*, *Sinh-Sinh* and *Exp-Sinh* quadrature speed and accuracy.  Included is the source code in C and BASIC (e.g. for calculators) to demonstrate wider applicability.

The two implementations `qthsh` *"cutiesh"* (Tanh-Sinh) and `quad` (Tanh-Sinh, Sinh-Sinh, Exp-Sinh) were extensively tested with 1084 integrals.

*Tanh-Sinh quadrature* is a method for numerical integration introduced by Hidetoshi Takahashi and Masatake Mori.  The method uses the tanh and sinh hyperbolic functions in a change of variable to transform the (−1,+1) open interval of the integral to an open interval on the entire real line (−∞,+∞).  Singularities at one or both endpoints of the (−1,+1) interval are mapped to the (−∞,+∞) endpoints of the transformed interval, forcing the endpoint singularities to vanish.  This makes the method quite insensitive to endpoint behavior, resulting in a significant enhancement of the accuracy of the numerical integration procedure compared to quadrature formulas that are based on the *trapezoidal* or *midpoint* rules with equidistant grids.  In most cases, the transformed integrand displays a rapid roll-off (decay) at a *double exponential* rate, enabling the numerical integrator to quickly achieve convergence.  This method is therefore also known as the *Double Exponential* (DE) formula.  A modification of the *Tanh-Sinh* formula was introduced by Krzysztof Michalski and Juan Mosig.  This modification simplifies the formulas for the abscissas and weights.  This modification requires fewer arithmetic operations to speed up numerical integration.

For more details on these methods and the optimizations applied, please see [Improving the Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas](https://www.genivia.com/files/qthsh.pdf).

## Examples

Level max `n=6` is recommended for IEEE 754 double floating point precision and a relative error tolerance `eps=1e-9` for maximum precision of the result with 9 digits or more:

    double f(double x) { return acos(x); }
    printf("%.15g\n", quad(f, 0, 1, 6, 1e-9, NULL));

Displays `1`, the exact integral with 15 digits precise.

    double f(double x) { return exp(-x/5); }
    printf("%.15g\n", quad(f, 0, INFINITY, 6, 1e-9, NULL));

Displays `5`, the exact integral with 15 digits precise.

    double f(double x) { return pow(cosh(x),-2); }
    printf("%.15g\n", quad(f, -INFINITY, INFINITY, 6, 1e-9, NULL));

Displays `2`, the exact integral with 15 digits precise.

#pragma once
#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>
#include <vector>

#include "RootForgeCore.h"


using namespace std;

#ifndef ROOTFORGE_POLYNOMIAL_H
#define ROOTFORGE_POLYNOMIAL_H

#endif //ROOTFORGE_POLYNOMIAL_H
class Polynomial
{
public:
    RootForge::Function f;
    vector<double> coeffs;
public:
    explicit Polynomial(const vector<double>& coefficients);
    int degree() const;
    RootForge::Function return_func() const;
    double evaluate(const double x) const;
    complex<double> evaluate(const complex<double>& x) const;
    static Polynomial add(const Polynomial p1, const Polynomial p2);
    static Polynomial negate(const Polynomial p1, Polynomial p2);
    static Polynomial ScalarProd(Polynomial p, const double x);
    static vector<Polynomial> divide(Polynomial p1, const Polynomial p2);
    static vector<double> multiply(Polynomial p1, Polynomial p2);
    static vector<double> NewtonHorner(Polynomial p);
    double NumericalDiff(const double x) const;
    static Polynomial Differentiate(Polynomial& p);
    static Polynomial Integrate(Polynomial& p);
    double NumericalInt(const double x_start,const double x_end) const;
    static void quadratic_formula(double a, double b, double c);
};

#pragma once

#include <cmath>
#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include "RootForgeCore.h"
#include "RootForge.h"
using namespace std;

#define PI 3.14159265359

#ifndef ROOTFORGE_MATHEMATICALFUNCTIONS_H
#define ROOTFORGE_MATHEMATICALFUNCTIONS_H

#endif //ROOTFORGE_MATHEMATICALFUNCTIONS_H

namespace RootForge::MathematicalFunctions
{
    int sgn(const double x);
    double FractionalPart(const double x);
    double power(const double a, const double n);
    bool is_integer(const double num);
    vector<int> sieve_of_eratosthenes(const int n);
    long double factorial(const long n);
    vector<complex<double>> nthRoot(int n, bool invert = false);
    double Gamma(const double n);
    double SpougeGamma(const double n);
    double LogGamma(const double n);
    long BinomialCoeff(const long n, const long k);
    double erf(const double x);
    double erfc(const double x);
    double taylorExp(const double a, const double x);
    double taylorSin(const double a, const double x);
    double taylorCos(const double a, const double x);
    double taylorLn(double x);
    string simplify_sqrt(int n);
    bool isPrime(int n);
    vector<complex<double>> FFT(vector<double> x_axis);
};

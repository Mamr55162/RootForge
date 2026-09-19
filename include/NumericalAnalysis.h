#pragma once
#include <iostream>
#include <cmath>
#include <vector>

#include "RootForgeCore.h"

using namespace std;

#ifndef ROOTFORGE_NUMERICALANALYSIS_H
#define ROOTFORGE_NUMERICALANALYSIS_H

#endif //ROOTFORGE_NUMERICALANALYSIS_H

namespace RootForge::Differentiation
{
    double forwardDifference(const Function& f, const double x, const double h);
    double backwardDifference(const Function& f, const double x, const double h);
    double centeredDifference(const Function& f, const double x, const double h);
    double D2centeredDifference(const Function& f, const double x, const double h);
    double RichardsonExtrapolationDiff(const Function& f, const double x, double h);
}
namespace RootForge::Integration
{
    double trapezoidal(const Function& f, const double start, const double end, const int n);
    double midpoint(const Function& f, const double start, const double end, const int n);
    double simpsons_13(const Function& f, const double start, const double end, const int n);
    double RichardsonExtrapolationInt(const Function& f, const double start, const double end, int n);
    double GaussianQuadrature(const Function& f, const double a, const double b);
}
namespace RootForge::Root
{
    double NewtonMethod(const double num, const double n);
    double Bisection(const Function& f, double start, double end);
    double NewtonRaphson(const Function& f);
    double Secant(const Function& f);
    double RegulaFalsi(const Function& f);
}

namespace RootForge::Interpolation
{
    double LinearInterpolation(double x1, double y1, double x2, double y2, double x_i);
    double NumericalPolyInterpolation(const vector<double>& x_axis, const vector<double>& y_axis, double x);
    vector<double> LagrangePolyInterpolation(const vector<double>& x_axis, const vector<double>& y_axis);
}

namespace RootForge::ODE
{
    double EulerMethod(const DFunction& f, const double x0, const double y0, const double h, const double target);
    double RK4(const DFunction& f, const double x0, const double y0, const double h, const double target);
    double Midpoint(const DFunction& f, const double x0, const double y0, const double h, const double target);
    double RK2(const DFunction& f, const double x0, const double y0, const double h, const double target);
    vector<double> SystemRK4(const vector<NFunction>& functions,const double x0, vector<double> y, const double h, const double target);
}

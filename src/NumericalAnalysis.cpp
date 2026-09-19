#include "NumericalAnalysis.h"
#include "Polynomial.h"
#include <stdexcept>

#include "MathematicalFunctions.h"

namespace RootForge::Differentiation
{
    double forwardDifference(const Function& f, const double x, const double h)
    {
        if (h == 0.0) return 0.0;
        return (f(x + h) - f(x)) / h;
    }
    double backwardDifference(const Function& f, const double x, const double h)
    {
        if (h == 0.0) return 0.0;
        return (f(x) - f(x - h)) / h;
    }
    double centeredDifference(const Function& f, const double x, const double h)
    {
        if (h == 0.0) return 0.0;
        return (f(x + h) - f(x - h)) / (2*h);
    }
    double D2centeredDifference(const Function& f, const double x, const double h)
    {
        if (h == 0.0) return 0.0;
        return (f(x + h) - 2*f(x) + f(x - h)) / (h*h);
    }
    double RichardsonExtrapolationDiff(const Function& f, const double x, double h)
    {
        const double Q_1 = centeredDifference(f, x, h);
        const double Q_2 = centeredDifference(f, x, h/2);
        return ( (4.0 * Q_2) - Q_1) / 3.0;
    }
}
namespace RootForge::Integration
{
    double trapezoidal(const Function& f, const double start, const double end, const int n)
    {
        double delta_x = (end - start) / n;
        double mid = 0;
        for (int i = 1; i <= n-1; i++)
            mid += f(i);
        return (delta_x / 2) * (f(start) + 2*mid + f(end));
    }
    double midpoint(const Function& f, const double start, const double end, const int n)
    {
        double delta_x = (end - start) / n;
        double sum = 0.0;
        for (int i = 1; i <= n; i++)
        {
            double m = start + (i - 0.5) * delta_x;
            sum += f(m);
        }
        return delta_x * sum;
    }
    double simpsons_13(const Function& f, const double start, const double end, const int n)
    {
        if (n % 2 != 0) return 0.0;
        double delta_x = (end - start) / n;
        double sum = 0.0;
        for (int i = 1; i <= n - 1; i++)
        {
            double xi = start + i*delta_x;
            if (i % 2 != 0)
                sum += 4*f(xi);
            else
                sum += 2*f(xi);
        }
        return (delta_x / 3) * (f(start) + sum + f(end));
    }
    double RichardsonExtrapolationInt(const Function& f, const double start, const double end, int n)
    {
        const double S_1 = simpsons_13(f, start, end, n);
        const double S_2 = simpsons_13(f, start, end, 2*n);
        return ( (16.0/15.0) * S_2) - ( (1.0/15.0) * S_1);
    }
    double GaussianQuadrature(const Function& f, const double a, const double b)
    {
        double front_term = (b - a) / 2.0;
        vector<double> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        vector<double> nodes = {-sqrt(0.6), 0, sqrt(0.6)};
        int n = 3;
        double result = 0.0;
        for (int i = 1; i <= n; i++)
        {
            result += weights[i - 1] * f( ((b-a) / 2.0) * nodes[i - 1] + ( (b + a) / 2.0) );
        }
        return (front_term * result);
    }
}
namespace RootForge::Root
{
    double NewtonMethod(const double num, const double n)
    {
        NumericalResult result;
        result.value = 0.5 * (num + 1);
        double prevValue;
        result.iterations = 0;
        do
        {
            prevValue = result.value;
            result.value = (1.0 / n) * ((n - 1) * prevValue + num / pow(prevValue, n - 1));
            result.iterations++;
        } while (fabs(prevValue - result.value) > 1E-9);
        const double actualRoot = pow(num, 1.0 / n);
        result.absError = fabs(actualRoot - result.value);
        result.relativeError = result.absError / actualRoot;
        result.percentError = result.relativeError * 100;
        NumericalResult::ShowResult(result);
        return result.value;
    }
    double Bisection(const Function& f, double start, double end)
    {
        if (f(start) * f(end) > 0) return 0;
        double midpoint;
        while (fabs(end - start) > 1E-9)
        {
            midpoint = (start + end) / 2;
            if (f(midpoint) == 0) return midpoint;
            else if (f(midpoint) * f(start) < 0)
            {
                end = midpoint;
            }
            else if (f(midpoint) * f(end) < 0)
            {
                start = midpoint;
            }
        }
        return midpoint;
    }
    double NewtonRaphson(const Function& f)
    {
        double x = 1.0;
        double guess = x - (f(x) / Differentiation::centeredDifference(f, x, 1E-5));
        while (fabs(guess - x) > 1E-9)
        {
            guess = x;
            x = x - (f(x) / Differentiation::centeredDifference(f, x, 1E-5));
        }
        return x;
    }
    double Secant(const Function& f)
    {
        double x0;
        double x1;
        for (int i = -100; i < 100; i++)
        {
            if (f(i) * f(i - 1) < 0)
            {
                x1 = i;
                x0 = i - 1;
                break;
            }
        }
        double x = 1;
        double guess = x1 - f(x1) * ((x1 - x0)/(f(x1) - f(x0)));
        while (fabs(guess - x) > 1E-9)
        {
            guess = x;
            x = x1 - f(x1) * ((x1 - x0)/(f(x1) - f(x0)));
            x0 = x1;
            x1 = x;
        }
        return x;
    }
    double RegulaFalsi(const Function& f)
    {
        double a = 0;
        double b = 0;
        for (int i = -100; i < 100; i++)
        {
            if (f(i) * f(i - 1) < 0)
            {
                b = i;
                a = i - 1;
                break;
            }
        }
        double c;
        while (fabs(f(c)) > 1E-9)
        {
            c = (a * f(b) - b * f(a)) / (f(b) - f(a));
            if (f(c) * f(a) < 0)
            {
                b = c;
            }
            else if (f(c) * f(b) < 0)
            {
                a = c;
            }
        }
        return c;
    }
}

namespace RootForge::Interpolation
{
    double LinearInterpolation(const double x1, const double y1, const double x2, const double y2, const double x_i)
    {
        if (x2 == x1) return NAN;
        return y1 + ( (x_i - x1) * (y2 - y1)/ (x2 - x1));
    }
    double NumericalPolyInterpolation(const vector<double>& x_axis, const vector<double>& y_axis, const double x)
    {
        if (x_axis.size() != y_axis.size()) return {0.0};
        double L = 1;
        double result = 0.0;
        const int n = x_axis.size();
        for (int i = 0; i < n; i++)
        {
            L = 1;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    L *= ((x - x_axis[j]) / (x_axis[i] - x_axis[j]));
            }
            result += (y_axis[i] * L);
        }
        return result;
    }
    vector<double> LagrangePolyInterpolation(const vector<double>& x_axis, const vector<double>& y_axis)
    {
        if (x_axis.size() != y_axis.size()) return {0.0};
        int degree = x_axis.size() - 1;
        vector<Polynomial> final;
        for (int i = 0; i <= degree; i++)
        {
            vector<Polynomial> poly;
            for (int j = 0; j <= degree; j++)
            {
                if (i != j)
                {
                    Polynomial p({(-x_axis[j] / (x_axis[i] - x_axis[j])), (1.0 / (x_axis[i] - x_axis[j]))});
                    poly.push_back(p);
                }
            }
            vector<double> coefficients = poly[0].coeffs;
            for (int k = 1; k < poly.size(); k++)
            {
                const Polynomial next(coefficients);
                coefficients = Polynomial::multiply(next, poly[k]);
            }
            for (int t = 0; t < coefficients.size(); t++)
            {
                coefficients[t] *= y_axis[i];
            }
            Polynomial current(coefficients);
            final.push_back(current);
        }
        Polynomial result = final[0];
        for (int i = 1; i < final.size(); i++)
        {
            result = Polynomial::add(result, final[i]);
        }
        return result.coeffs;
    }
}

namespace RootForge::ODE
{
    double EulerMethod(const DFunction& f, const double x0, const double y0, const double h, const double target)
    {
        double x = x0;
        double y = y0;
        double slope;
        while (x < target - 1e-9)
        {
            slope = f(x, y);
            x = x + h;
            y = y + h * slope;
        }
        return y;
    }
    double RK4(const DFunction& f, const double x0, const double y0, const double h, const double target)
    {
        double x = x0;
        double y = y0;
        double k1;
        double k2;
        double k3;
        double k4;
        double current_h = h;
        while (x < target - 1E-9)
        {
            if (x + current_h > target)
                current_h = target - x;

            k1 = current_h * f(x, y);
            k2 = current_h * f(x + current_h/2, y +  k1/2);
            k3 = current_h * f(x + current_h/2, y +  k2/2);
            k4 = current_h * f(x + current_h, y + k3);
            x = x + current_h;
            y = y + (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);
        }
        return y;
    }
    double Midpoint(const DFunction& f, const double x0, const double y0, const double h, const double target)
    {
        double x = x0;
        double y = y0;
        while (x < target)
        {

            double current_h = (x + h > target) ? (target - x) : h;

            double x_mid = x + (current_h / 2);
            double y_mid = y + (current_h / 2) * f(x, y);

            y = y + current_h * f(x_mid, y_mid);
            x += current_h;
        }

        return y;
    }
    double RK2(const DFunction& f, const double x0, const double y0, const double h, const double target)
    {
        double x = x0;
        double y = y0;

        double k1, k2;
        while (x < target - 1E-9)
        {
            double current_h = (x + h > target) ? (target - x) : h;
            k1 = f(x, y);
            k2 = f(x + current_h, y + (current_h * k1));
            y += (current_h / 2.0) * (k1 + k2);
            x += current_h;
        }
        return y;
    }
    vector<double> SystemRK4(const vector<NFunction>& functions, const double x0, vector<double> y, const double h, const double target)
    {
        double x = x0;
        double current_h = h;
        const size_t n = y.size();
        vector<double> y_temp(n);
        vector<double> k1(n), k2(n), k3(n), k4(n);
        vector<double> midpoint_1(n), midpoint_2(n), midpoint_3(n);
        while (x < target - 1E-9)
        {
            if (x + current_h > target)
                current_h = target - x;


            for (int i = 0; i < n; i++)
            {
                k1[i] = (functions[i](x, y));
                midpoint_1[i] = y[i] + (current_h / 2.0) * k1[i];
            }

            for (int i = 0; i < n; i++)
            {
                k2[i] = (functions[i](x + (current_h / 2.0), midpoint_1));
                midpoint_2[i] = y[i] + (current_h / 2.0) * k2[i];
            }

            for (int i = 0; i < n; i++)
            {
                k3[i] = (functions[i](x + (current_h / 2.0), midpoint_2));
                midpoint_3[i] = y[i] + current_h * k3[i];
            }

            for (int i = 0; i < n; i++)
            {
                k4[i] = (functions[i](x + (current_h), midpoint_3));
            }
            for (int i = 0; i < n; i++)
            {
                y[i] += (current_h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            x += current_h;
        }
        return y;
    }
}

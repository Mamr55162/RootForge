#include "Polynomial.h"
#include "MathematicalFunctions.h"
#include "NumericalAnalysis.h"
#include <stdexcept>

Polynomial::Polynomial(const vector<double>& coefficients) : coeffs(coefficients)
{
    f = [coefficients](const double x)
    {
        double result = coefficients.back();
        for (int i = coefficients.size() - 2; i >= 0; i--)
        {
            result = result * x + coefficients[i];
        }
        return result;
    };
}
int Polynomial::degree() const
{
    return coeffs.size() - 1;
}
RootForge::Function Polynomial::return_func() const
{
    return f;
}
double Polynomial::evaluate(const double x) const
{
    if (coeffs.empty()) return 0.0;

    const vector<double> cof = coeffs;
    double term = cof[cof.size() - 1];
    for (int i = cof.size() - 2; i >= 0; i--)
    {
        term = (term * x) + cof[i];
    }
    return term;
}
complex<double> Polynomial::evaluate(const complex<double>& x) const
{
    if (coeffs.empty()) return complex<double>(0.0, 0.0);

    complex<double> term = coeffs.back();
    for (int i = coeffs.size() - 2; i >= 0; i--)
    {
        term = (term * x) + complex<double>(coeffs[i], 0.0);
    }
    return term;
}
Polynomial Polynomial::add(const Polynomial p1, const Polynomial p2)
{
    Polynomial p3({});
    size_t size = max(p1.coeffs.size(), p2.coeffs.size());
    p3.coeffs.resize(size, 0);
    for (int i = 0; i < p1.coeffs.size(); i++)
        p3.coeffs[i] += p1.coeffs[i];
    for (int i = 0; i < p2.coeffs.size(); i++)
        p3.coeffs[i] += p2.coeffs[i];
    return p3;
}
Polynomial Polynomial::negate(const Polynomial p1, Polynomial p2)
{
    for (int i = 0; i < p2.coeffs.size(); i++)
        p2.coeffs[i] *= -1;
    return add(p1, p2);
}
Polynomial Polynomial::ScalarProd(Polynomial p, const double x)
{
    for (int i = 0; i < p.coeffs.size(); i++)
        p.coeffs[i] *= x;
    return p;
}
vector<Polynomial> Polynomial::divide(Polynomial p1, const Polynomial p2)
{
    vector<double> quotient_coeff;
    Polynomial p3({});
    p3.coeffs.resize(p1.coeffs.size(), 0);
    int deg1 = p1.degree();
    const int deg2 = p2.degree();
    if (deg1 < deg2)
        return {Polynomial({0.0}), p1};
    int deg3 = p3.degree();
    while (deg1 >= deg2)
    {
        fill(p3.coeffs.begin(), p3.coeffs.end(),0);
        double term = p1.coeffs[deg1] / p2.coeffs[deg2];
        quotient_coeff.push_back(term);
        for (int i = 0; i < p2.coeffs.size(); i++)
        {
            p3.coeffs[deg3 - i] = term * p2.coeffs[deg2 - i];
        }
        p1 = negate(p1, p3);
        deg1--;
        deg3--;
    }
    reverse(quotient_coeff.begin(), quotient_coeff.end());
    Polynomial result(quotient_coeff);
    return {result, p1};
}
vector<double> Polynomial::multiply(Polynomial p1, Polynomial p2)
{
    int n = RootForge::MathematicalFunctions::power(p1.degree() + p2.degree(), 2);
    vector<double> A = p1.coeffs;
    vector<double> B = p2.coeffs;
    A.resize(n, 0.0);
    B.resize(n, 0.0);
    vector<complex<double>> unity_roots = RootForge::MathematicalFunctions::nthRoot(n);
    complex<double> A_star;
    complex<double> B_star;
    vector<complex<double>> C(n);
    for (int i = 0; i < n; i++)
    {
        A_star = p1.evaluate(unity_roots[i]);
        B_star = p2.evaluate(unity_roots[i]);
        C[i] = (A_star * B_star);
    }
    vector<complex<double>> result;
    vector<complex<double>> reverse_roots = RootForge::MathematicalFunctions::nthRoot(n, true);
    complex<double> term = 0.0;
    for (int i = 0; i < n; i++)
    {
        term = 0.0;
        for (int j = 0; j < n; j++)
        {
            term += C[j] * reverse_roots[(i * j) % n];
        }
        result.push_back((1.0/n) * term);
    }
    vector<double> final_result(result.size());
    for (int i = 0; i < result.size(); i++)
        final_result[i] = (result[i].real());
    return final_result;
}
vector<double> Polynomial::NewtonHorner(Polynomial p)
{
    vector<double> result;
    int degree = p.degree();
    const double EPSILON = 1e-9;
    while (degree > 0)
    {
        RootForge::Function current_polynomial = p.return_func();
        double root = RootForge::Root::NewtonRaphson(current_polynomial);
        result.push_back(root);
        Polynomial current_root({-root, 1});
        vector<Polynomial> res = divide(p, current_root);
        p = res[0];
        degree--;
        if (degree == 1 && fabs(p.coeffs[0]) < EPSILON) break;
    }
    return result;
}
double Polynomial::NumericalDiff(const double x) const
{
    return RootForge::Differentiation::centeredDifference(f, x, 1E-5);
}
Polynomial Polynomial::Differentiate(Polynomial& p)
{
    if (p.coeffs.size() <= 1) {
        return Polynomial({0.0});
    }
    vector<double> derived_coeffs;
    for (size_t i = 1; i < p.coeffs.size(); i++)
    {
        derived_coeffs.push_back(p.coeffs[i] * i);
    }
    p.coeffs = derived_coeffs;
    return p;
}
Polynomial Polynomial::Integrate(Polynomial& p)
{
    if (p.coeffs.size() == 0) return Polynomial({0.0});
    vector<double> integrated_coeffs;
    integrated_coeffs.push_back(0.0);
    for (int i = 0; i < p.coeffs.size(); i++)
    {
        integrated_coeffs.push_back(p.coeffs[i] / (i + 1));
    }
    p.coeffs = integrated_coeffs;
    return p;
}
double Polynomial::NumericalInt(const double x_start,const double x_end) const
{
    return RootForge::Integration::simpsons_13(f, x_start, x_end, 100);
}
void Polynomial::quadratic_formula(double a, double b, double c)
{
    if (a == 0)
    {
        cerr << "a cannot equal 0\n";
        return;
    }
    double discrimant = pow(b, 2) - 4*a*c;
    if (discrimant > 0)
    {
        double sol1 = (-b - sqrt(discrimant)) / (2*a);
        double sol2 = (-b + sqrt(discrimant)) / (2*a);
        cout << "Root 1 = " << sol1 << endl;
        cout << "Root 2 = " << sol2 << endl;
    }
    else if (discrimant == 0)
    {
        double sol = -b / (2*a);
        cout << "Root = " << sol << endl;
    }
    else
    {
        double sol_real = -b / (2*a);
        double sol_imag = (-discrimant) / (4*pow(a,2));
        string imag = RootForge::MathematicalFunctions::simplify_sqrt(sol_imag);
        cout << "Root1 = " << sol_real << " + " << imag << "i" << endl;
        cout << "Root2 = " << sol_real << " - " << imag << "i";
    }
}

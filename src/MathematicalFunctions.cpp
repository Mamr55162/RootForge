#include "MathematicalFunctions.h"
#include "NumericalAnalysis.h"
#include <stdexcept>

namespace RootForge::MathematicalFunctions
{
    int sgn(const double x)
    {
        if (x < 0) return -1;
        else if (x == 0) return 0;
        else return 1;
    }
    double FractionalPart(const double x)
    {
        return x - floor(x);
    }
    double power(const double a, const double n)
    {
        bool fractional = (fabs(n - floor(n)) != 0);
        if (fractional) return exp(n * log(a));
        long long exp = static_cast<int> (n);
        if (a == 0.0 && exp < 0) {
            return INFINITY;
        }
        double result = 1.0;
        double base = a;
        if (exp < 0)
        {
            base = 1.0 / base;
            exp = -exp;
        }
        while (exp > 0)
        {
            if (exp % 2 != 0)
                result *= base;
            base *= base;
            exp /= 2;
        }
        return result;
    }
    bool is_integer(const double num)
    {
        return floor(num) == num;
    }
    vector<int> sieve_of_eratosthenes(const int n)
    {
        bool isPrime[n + 1];
        for (int i = 2; i < n + 1; i++)
            isPrime[i] = true;
        isPrime[0] = false;
        isPrime[1] = isPrime[0];
        for (int i = 2; i * i < n; i++)
        {
            for (int j = i * i; j < n + 1; j++)
            {
                if (j % i == 0 && isPrime[j] != false)
                    isPrime[j] = false;
            }
        }
        vector<int> result;
        for (int i = 0; i < n + 1; i++)
        {
            if (isPrime[i])
                result.push_back(i);
        }
        return result;
    }
    long double factorial(const long n)
    {
        long double result = 1;
        const vector<int> primes = sieve_of_eratosthenes(n);
        vector<int> exponents;
        for (int i = 0; i < primes.size(); i++)
        {
            int current_exponent = 0;
            int k = 1;
            while (floor(n / (power(primes[i], k))) > 0)
            {
                current_exponent += floor(n / power(primes[i], k));
                k++;
            }
            exponents.push_back(current_exponent);
        }
        for (int i = 0; i < primes.size(); i++)
        {
            result *= power(primes[i], exponents[i]);
        }
        return result;
    }
    vector<complex<double>> nthRoot(int n, bool invert)
    {
        vector<complex<double>> roots;
        roots.reserve(n);
        complex<double> z;
        double sign = invert ? 1.0 : -1.0;
        for (int i = 0; i < n; i++)
        {
            double theta = (sign * 2.0 * PI * i) / n;
            z = polar(1.0, theta);
            roots.push_back(z);
        }
        return roots;
    }
    double Gamma(const double n)
    {
        if (n <= 0 && is_integer(n)) return NAN;
        if (n > 0 && is_integer(n)) return factorial(static_cast<int>(n) - 1);
        if (n == 0.5) return sqrt(PI);
        bool is_half_integer = is_integer(n - 0.5);
        if (!is_half_integer) return NAN;
        if (n > 0) return (n - 1) * Gamma(n - 1);
        return  Gamma(n + 1) / n;
    }
    double SpougeGamma(const double n)
    {
        double result = 0;
        double sum = 0;
        const double z = n - 1;
        const int a = 3;
        double c_k;
        const double c0 = sqrt(2*PI);
        const double front = power(z + a, z + 0.5) * exp(-(z + a));
        for (int k = 1; k < a; k++)
        {
            c_k = (power(-1, k-1) / factorial(k-1)) * (power(a-k, k - 0.5)) * exp(a - k);
            sum += (c_k) / (z + k);
        }
        result = front * (c0 + sum);
        return result;
    }
    double LogGamma(const double n)
    {
        return log(Gamma(n));
    }
    long BinomialCoeff(const long n, const long k)
    {
        if (n < k) return 0;
        if (n * k < 0) return 0;
        return (factorial(n)) / (factorial(k) * factorial(n - k));
    }
    double erf(const double x)
    {
        if (x > 4.0) return 1.0;
        if (x < -4.0) return -1.0;
        return (2 / sqrt(PI)) * (Integration::simpsons_13([](const double t){return exp(-t*t);}, 0, x, 100));
    }
    double erfc(const double x)
    {
        const double result = 1.0 - erf(x);
        return (result < 0.0) ? 0.0 : result;
    }
//============== Exponential (e^ax) Taylor Expansion ==============
    double taylorExp(const double a, const double x)
    {
        NumericalResult result;
        if (a == 0 || x == 0) return 1;
        result.value = 1;
        int n = 1;
        double term = 1;
        result.iterations = 0;
        while (fabs(term) > 1E-18)
        {
            term *= (a*x) / n;
            result.value += term;
            n++;
            result.iterations++;
        }
        result.absError = fabs(exp(a*x) - result.value);
        result.relativeError = result.absError / exp(a*x);
        result.percentError = 100 * result.relativeError;
        NumericalResult::ShowResult(result);
        return result.value;
    }
//============== Sine (sin(ax)) Taylor Expansion ==============
    double taylorSin(const double a, const double x)
    {
        NumericalResult result;
        if (a == 0 || x == 0) return 0;
        result.value = 0;
        int n = 1;
        double term = a * x;
        result.iterations = 0;
        while (fabs(term) > 1E-9)
        {
            result.value += term;
            term = -term * ((a*x) * (a*x)) / ((2 * n) * (2 * n + 1));
            n++;
            result.iterations++;
        }
        result.absError = fabs(sin(a*x) - result.value);
        result.relativeError = result.absError / sin(a*x);
        result.percentError = 100 * result.relativeError;
        NumericalResult::ShowResult(result);
        return result.value;
    }
//============== Cosine (cos(ax)) Taylor Expansion ==============
    double taylorCos(const double a, const double x)
    {
        NumericalResult result;
        if (a == 0 || x == 0) return 1;
        result.value = 0;
        int n = 0;
        double term = 1;
        result.iterations = 0;
        while (fabs(term) > 1E-9)
        {
            result.value += term;
            term = ( (-term) * ((a*x) * (a*x)) ) / ((2*n + 2) * (2*n + 1));
            n++;
            result.iterations++;
        }
        result.absError = fabs(cos(a*x) - result.value);
        result.relativeError = result.absError / cos(a*x);
        result.percentError = 100 * result.relativeError;
        NumericalResult::ShowResult(result);
        return result.value;
    }
//============== Natural Logarithm (ln(x)) Taylor Expansion ==============
    double taylorLn(double x)
    {
        NumericalResult result;
        if (x == 0)
        {
            cerr << "Undefined\n";
            return -1;
        }
        if (x == 1) return 0;
        result.value = 0.0;
        int n = 1;
        double term = (x - 1) / (x + 1);
        result.iterations = 0;
        while (fabs(term) > 1E-9)
        {
            result.value += term;
            term = ( term *  (2*n + 1) * ((x-1)/(x+1)) * ((x-1)/(x+1)) ) / (2*n + 3) ;
            n++;
            result.iterations++;
        }
        result.value *= 2;
        result.absError = fabs(log(x) - result.value);
        result.relativeError = result.absError / log(x);
        result.percentError = 100 * result.relativeError;
        NumericalResult::ShowResult(result);
        return result.value;
    }
//============== Simplify square root ==============
    string simplify_sqrt(int n)
    {
        if (n < 0) return "Invalid input";
        if (n == 0) return"0";
        if (static_cast<int> (sqrt(n)) == sqrt(n))
        {
            return to_string(static_cast<int> (sqrt(n)));
        }
        string result;

        for (int i = sqrt(n); i > 1; i--)
        {
            if ((n) % (i * i) == 0)
            {

                int outer = i;
                int inner = n / (i * i);
                result += to_string(outer);
                result += "\xE2\x88\x9A";
                result += to_string(inner);
                return result;

            }
        }
        return "\xE2\x88\x9A" + std::to_string(n);
    }

    bool isPrime(int n)
    {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        for (int i = 5; i * i <= n; i += 6)
        {
            if (n % i == 0 || n % (i + 2) == 0) return false;
        }
        return true;
    }
}

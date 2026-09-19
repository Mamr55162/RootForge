#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>

#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "RootForgeCore.h"
#include "MathematicalFunctions.h"
#include "NumericalAnalysis.h"
#include "Polynomial.h"

using namespace std;

// --------------------------------------------------------------------------
// Small helpers
// --------------------------------------------------------------------------

static void glfw_error_callback(int error, const char* description) {
    std::cerr << "GLFW Error " << error << ": " << description << std::endl;
}

// Parses "1, 2.5, -3" -> {1.0, 2.5, -3.0}. Silently skips anything that
// doesn't parse as a number.
static vector<double> ParseDoubles(const string& s)
{
    vector<double> out;
    stringstream ss(s);
    string item;
    while (getline(ss, item, ','))
    {
        try { out.push_back(stod(item)); } catch (...) { /* skip bad token */ }
    }
    return out;
}

static string FormatVec(const vector<double>& v)
{
    ostringstream oss;
    oss << "[ ";
    for (size_t i = 0; i < v.size(); i++)
    {
        oss << v[i];
        if (i + 1 < v.size()) oss << ", ";
    }
    oss << " ]";
    return oss.str();
}

static string FormatComplexVec(const vector<complex<double>>& v)
{
    ostringstream oss;
    oss << fixed << setprecision(4);
    for (size_t i = 0; i < v.size(); i++)
    {
        oss << v[i].real();
        if (v[i].imag() >= 0) oss << " + " << v[i].imag() << "i";
        else oss << " - " << fabs(v[i].imag()) << "i";
        if (i + 1 < v.size()) oss << "\n";
    }
    return oss.str();
}

// --------------------------------------------------------------------------
// Predefined test functions, since there is no expression parser in the lib
// --------------------------------------------------------------------------

static const char* kTestFuncNames[] = {
    "f(x) = x^2 - 4",
    "f(x) = sin(x)",
    "f(x) = cos(x)",
    "f(x) = e^x - 2",
    "f(x) = x^3 - x - 2",
    "f(x) = ln(x) - 1"
};

static RootForge::Function GetTestFunction(int idx)
{
    switch (idx)
    {
        case 0: return [](double x){ return x * x - 4.0; };
        case 1: return [](double x){ return sin(x); };
        case 2: return [](double x){ return cos(x); };
        case 3: return [](double x){ return exp(x) - 2.0; };
        case 4: return [](double x){ return x * x * x - x - 2.0; };
        case 5: return [](double x){ return log(x) - 1.0; };
        default: return [](double x){ return x; };
    }
}

static const char* kODEFuncNames[] = {
    "dy/dx = y",
    "dy/dx = x + y",
    "dy/dx = x^2 - y",
    "dy/dx = cos(x)"
};

static RootForge::DFunction GetODEFunction(int idx)
{
    switch (idx)
    {
        case 0: return [](double x, double y){ (void)x; return y; };
        case 1: return [](double x, double y){ return x + y; };
        case 2: return [](double x, double y){ return x * x - y; };
        case 3: return [](double x, double y){ (void)y; return cos(x); };
        default: return [](double x, double y){ (void)x; (void)y; return 0.0; };
    }
}

// --------------------------------------------------------------------------
// Library / function menu definitions
// --------------------------------------------------------------------------

static const char* kMathFuncs[] = {
    "sgn(x)", "FractionalPart(x)", "power(a, n)", "is_integer(x)",
    "sieve_of_eratosthenes(n)", "factorial(n)", "nthRoot(n, invert)",
    "Gamma(n)", "SpougeGamma(n)", "LogGamma(n)", "BinomialCoeff(n, k)",
    "erf(x)", "erfc(x)", "taylorExp(a, x)", "taylorSin(a, x)",
    "taylorCos(a, x)", "taylorLn(x)", "simplify_sqrt(n)", "isPrime(n)"
};

static const char* kDiffFuncs[] = {
    "forwardDifference(f, x, h)", "backwardDifference(f, x, h)",
    "centeredDifference(f, x, h)", "D2centeredDifference(f, x, h)",
    "RichardsonExtrapolationDiff(f, x, h)"
};

static const char* kIntFuncs[] = {
    "trapezoidal(f, a, b, n)", "midpoint(f, a, b, n)", "simpsons_13(f, a, b, n)",
    "RichardsonExtrapolationInt(f, a, b, n)", "GaussianQuadrature(f, a, b)"
};

static const char* kRootFuncs[] = {
    "NewtonMethod - nth root of a number", "Bisection(f, start, end)",
    "NewtonRaphson(f)", "Secant(f)", "RegulaFalsi(f)"
};

static const char* kInterpFuncs[] = {
    "LinearInterpolation(x1,y1,x2,y2,xi)", "NumericalPolyInterpolation(x,y,xi)",
    "LagrangePolyInterpolation(x,y)"
};

static const char* kODEFuncsList[] = {
    "EulerMethod", "RK4", "Midpoint", "RK2", "SystemRK4 (2-eq SHM demo)"
};

static const char* kPolyFuncs[] = {
    "Evaluate p(x)", "Add  p1 + p2", "Subtract  p1 - p2", "ScalarProd  k * p",
    "Divide  p1 / p2", "Multiply  p1 * p2", "NewtonHorner (find all real roots)",
    "Differentiate", "Integrate", "NumericalDiff(x)", "NumericalInt(x1, x2)",
    "Quadratic Formula (a, b, c)"
};

struct LibraryInfo { const char* name; const char* const* funcs; int count; };

static const LibraryInfo kLibraries[] = {
    { "Mathematical Functions", kMathFuncs,     IM_ARRAYSIZE(kMathFuncs) },
    { "Differentiation",        kDiffFuncs,     IM_ARRAYSIZE(kDiffFuncs) },
    { "Integration",            kIntFuncs,      IM_ARRAYSIZE(kIntFuncs) },
    { "Root Finding",           kRootFuncs,     IM_ARRAYSIZE(kRootFuncs) },
    { "Interpolation",          kInterpFuncs,   IM_ARRAYSIZE(kInterpFuncs) },
    { "ODE Solvers",            kODEFuncsList,  IM_ARRAYSIZE(kODEFuncsList) },
    { "Polynomial",             kPolyFuncs,     IM_ARRAYSIZE(kPolyFuncs) }
};
static const int kLibraryCount = IM_ARRAYSIZE(kLibraries);

int main() {
    // 1. Initialize GLFW and create a window
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return 1;

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(1000, 780, "RootForge Library Explorer", nullptr, nullptr);
    if (window == nullptr) return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // 2. Initialize Dear ImGui Context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    ImGui::StyleColorsDark(); // Apply dark theme

    // 3. Initialize Platform/Renderer Backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // ---------------- Persistent UI / calculation state ----------------
    int currentLibrary  = 0;
    int currentFunction = 0;

    float valA = 1.0f, valB = 1.0f, valC = 1.0f, valD = 1.0f, valE = 0.0f;
    int   intN  = 10;
    bool  boolFlag  = false;
    int   testFuncIdx = 0;
    int   odeFuncIdx  = 0;

    char vecBufX[256]  = "1,2,3,4";
    char vecBufY[256]  = "1,4,9,16";
    char polyBufA[256] = "1,2,3";
    char polyBufB[256] = "1,1";

    string resultText = "Select a library and a function above, set the inputs, then press Calculate.";
    vector<string> history;

    // 4. The Main Application Loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- YOUR UI LOGIC STARTS HERE ---
        ImGui::SetNextWindowSize(ImVec2(960, 740), ImGuiCond_FirstUseEver);
        ImGui::Begin("RootForge Library Explorer");

        // --- Level 1: Library menu ---
        const char* libNames[kLibraryCount];
        for (int i = 0; i < kLibraryCount; i++) libNames[i] = kLibraries[i].name;

        if (ImGui::Combo("Library", &currentLibrary, libNames, kLibraryCount))
            currentFunction = 0; // reset the function choice whenever the library changes

        // --- Level 2: Function menu (depends on the chosen library) ---
        ImGui::Combo("Function", &currentFunction, kLibraries[currentLibrary].funcs,
                      kLibraries[currentLibrary].count);

        ImGui::Separator();
        ImGui::Spacing();

        // ------------------------------------------------------------------
        // Library 0: Mathematical Functions
        // ------------------------------------------------------------------
        if (currentLibrary == 0)
        {
            switch (currentFunction)
            {
                case 0: // sgn(x)
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        int r = RootForge::MathematicalFunctions::sgn(valA);
                        ostringstream oss; oss << "sgn(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 1: // FractionalPart
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::FractionalPart(valA);
                        ostringstream oss; oss << "FractionalPart(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 2: // power
                {
                    ImGui::InputFloat("a (base)", &valA);
                    ImGui::InputFloat("n (exponent)", &valB);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::power(valA, valB);
                        ostringstream oss; oss << "power(" << valA << ", " << valB << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 3: // is_integer
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        bool r = RootForge::MathematicalFunctions::is_integer(valA);
                        resultText = "is_integer(" + to_string(valA) + ") = " + (r ? "true" : "false");
                        history.push_back(resultText);
                    }
                    break;
                }
                case 4: // sieve_of_eratosthenes
                {
                    ImGui::InputInt("n", &intN);
                    if (ImGui::Button("Calculate"))
                    {
                        if (intN < 2) { resultText = "Enter n >= 2."; }
                        else
                        {
                            vector<int> primes = RootForge::MathematicalFunctions::sieve_of_eratosthenes(intN);
                            ostringstream oss; oss << "Primes up to " << intN << ": [ ";
                            for (size_t i = 0; i < primes.size(); i++) { oss << primes[i]; if (i + 1 < primes.size()) oss << ", "; }
                            oss << " ]";
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 5: // factorial
                {
                    ImGui::InputInt("n", &intN);
                    if (ImGui::Button("Calculate"))
                    {
                        long double r = RootForge::MathematicalFunctions::factorial(intN);
                        ostringstream oss; oss << intN << "! = " << (double)r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 6: // nthRoot
                {
                    ImGui::InputInt("n", &intN);
                    ImGui::Checkbox("invert", &boolFlag);
                    if (ImGui::Button("Calculate"))
                    {
                        if (intN <= 0) { resultText = "Enter n >= 1."; }
                        else
                        {
                            auto roots = RootForge::MathematicalFunctions::nthRoot(intN, boolFlag);
                            resultText = to_string(intN) + "-th roots of unity:\n" + FormatComplexVec(roots);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 7: // Gamma
                {
                    ImGui::InputFloat("n", &valA);
                    ImGui::TextDisabled("(Exact for integers and half-integers; NAN otherwise)");
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::Gamma(valA);
                        ostringstream oss; oss << "Gamma(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 8: // SpougeGamma
                {
                    ImGui::InputFloat("n", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::SpougeGamma(valA);
                        ostringstream oss; oss << "SpougeGamma(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 9: // LogGamma
                {
                    ImGui::InputFloat("n", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::LogGamma(valA);
                        ostringstream oss; oss << "LogGamma(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 10: // BinomialCoeff
                {
                    ImGui::InputInt("n", &intN);
                    static int kVal = 2;
                    ImGui::InputInt("k", &kVal);
                    if (ImGui::Button("Calculate"))
                    {
                        long r = RootForge::MathematicalFunctions::BinomialCoeff(intN, kVal);
                        ostringstream oss; oss << "C(" << intN << ", " << kVal << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 11: // erf
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::erf(valA);
                        ostringstream oss; oss << "erf(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 12: // erfc
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::erfc(valA);
                        ostringstream oss; oss << "erfc(" << valA << ") = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 13: // taylorExp
                {
                    ImGui::InputFloat("a", &valA);
                    ImGui::InputFloat("x", &valB);
                    ImGui::TextDisabled("(Full error report is also printed to the console)");
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::taylorExp(valA, valB);
                        ostringstream oss; oss << "taylorExp(a=" << valA << ", x=" << valB << ") ~= " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 14: // taylorSin
                {
                    ImGui::InputFloat("a", &valA);
                    ImGui::InputFloat("x", &valB);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::taylorSin(valA, valB);
                        ostringstream oss; oss << "taylorSin(a=" << valA << ", x=" << valB << ") ~= " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 15: // taylorCos
                {
                    ImGui::InputFloat("a", &valA);
                    ImGui::InputFloat("x", &valB);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::taylorCos(valA, valB);
                        ostringstream oss; oss << "taylorCos(a=" << valA << ", x=" << valB << ") ~= " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 16: // taylorLn
                {
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::MathematicalFunctions::taylorLn(valA);
                        ostringstream oss; oss << "taylorLn(" << valA << ") ~= " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 17: // simplify_sqrt
                {
                    ImGui::InputInt("n", &intN);
                    if (ImGui::Button("Calculate"))
                    {
                        string r = RootForge::MathematicalFunctions::simplify_sqrt(intN);
                        resultText = "simplify_sqrt(" + to_string(intN) + ") = " + r;
                        history.push_back(resultText);
                    }
                    break;
                }
                case 18: // isPrime
                {
                    ImGui::InputInt("n", &intN);
                    if (ImGui::Button("Calculate"))
                    {
                        bool r = RootForge::MathematicalFunctions::isPrime(intN);
                        resultText = "isPrime(" + to_string(intN) + ") = " + (r ? "true" : "false");
                        history.push_back(resultText);
                    }
                    break;
                }
            }
        }
        // ------------------------------------------------------------------
        // Library 1: Differentiation
        // ------------------------------------------------------------------
        else if (currentLibrary == 1)
        {
            ImGui::Combo("Test function f(x)", &testFuncIdx, kTestFuncNames, IM_ARRAYSIZE(kTestFuncNames));
            ImGui::InputFloat("x", &valA);
            ImGui::InputFloat("h (step size)", &valB, 0.0f, 0.0f, "%.6f");
            if (ImGui::Button("Calculate"))
            {
                RootForge::Function f = GetTestFunction(testFuncIdx);
                double r = 0.0; string name;
                switch (currentFunction)
                {
                    case 0: r = RootForge::Differentiation::forwardDifference(f, valA, valB); name = "forwardDifference"; break;
                    case 1: r = RootForge::Differentiation::backwardDifference(f, valA, valB); name = "backwardDifference"; break;
                    case 2: r = RootForge::Differentiation::centeredDifference(f, valA, valB); name = "centeredDifference"; break;
                    case 3: r = RootForge::Differentiation::D2centeredDifference(f, valA, valB); name = "D2centeredDifference"; break;
                    case 4: r = RootForge::Differentiation::RichardsonExtrapolationDiff(f, valA, valB); name = "RichardsonExtrapolationDiff"; break;
                }
                ostringstream oss;
                oss << name << " of " << kTestFuncNames[testFuncIdx] << " at x=" << valA << ", h=" << valB << " = " << r;
                resultText = oss.str(); history.push_back(resultText);
            }
        }
        // ------------------------------------------------------------------
        // Library 2: Integration
        // ------------------------------------------------------------------
        else if (currentLibrary == 2)
        {
            ImGui::Combo("Test function f(x)", &testFuncIdx, kTestFuncNames, IM_ARRAYSIZE(kTestFuncNames));
            ImGui::InputFloat("start", &valA);
            ImGui::InputFloat("end", &valB);
            if (currentFunction != 4) // GaussianQuadrature has no n
                ImGui::InputInt("n (subintervals)", &intN);
            if (currentFunction == 2 || currentFunction == 3)
                ImGui::TextDisabled("(simpsons_13 requires an even n)");

            if (ImGui::Button("Calculate"))
            {
                RootForge::Function f = GetTestFunction(testFuncIdx);
                if (currentFunction == 2 && intN % 2 != 0)
                {
                    resultText = "simpsons_13 requires an even n.";
                }
                else
                {
                    double r = 0.0; string name;
                    switch (currentFunction)
                    {
                        case 0: r = RootForge::Integration::trapezoidal(f, valA, valB, intN); name = "trapezoidal"; break;
                        case 1: r = RootForge::Integration::midpoint(f, valA, valB, intN); name = "midpoint"; break;
                        case 2: r = RootForge::Integration::simpsons_13(f, valA, valB, intN); name = "simpsons_13"; break;
                        case 3: r = RootForge::Integration::RichardsonExtrapolationInt(f, valA, valB, intN); name = "RichardsonExtrapolationInt"; break;
                        case 4: r = RootForge::Integration::GaussianQuadrature(f, valA, valB); name = "GaussianQuadrature"; break;
                    }
                    ostringstream oss;
                    oss << name << " of " << kTestFuncNames[testFuncIdx] << " on [" << valA << ", " << valB << "]";
                    if (currentFunction != 4) oss << " with n=" << intN;
                    oss << " = " << r;
                    resultText = oss.str();
                }
                history.push_back(resultText);
            }
        }
        // ------------------------------------------------------------------
        // Library 3: Root Finding
        // ------------------------------------------------------------------
        else if (currentLibrary == 3)
        {
            switch (currentFunction)
            {
                case 0: // NewtonMethod - nth root
                {
                    ImGui::InputFloat("number", &valA);
                    ImGui::InputFloat("root degree (n)", &valB);
                    ImGui::TextDisabled("(Full iteration report is also printed to the console)");
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::Root::NewtonMethod(valA, valB);
                        ostringstream oss; oss << valB << "-th root of " << valA << " = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 1: // Bisection
                {
                    ImGui::Combo("Test function f(x)", &testFuncIdx, kTestFuncNames, IM_ARRAYSIZE(kTestFuncNames));
                    ImGui::InputFloat("start", &valA);
                    ImGui::InputFloat("end", &valB);
                    if (ImGui::Button("Calculate"))
                    {
                        RootForge::Function f = GetTestFunction(testFuncIdx);
                        if (f(valA) * f(valB) > 0)
                        {
                            resultText = "No sign change on [" + to_string(valA) + ", " + to_string(valB) + "] - pick a different interval.";
                        }
                        else
                        {
                            double r = RootForge::Root::Bisection(f, valA, valB);
                            ostringstream oss; oss << "Root of " << kTestFuncNames[testFuncIdx] << " on [" << valA << ", " << valB << "] = " << r;
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 2: // NewtonRaphson
                case 3: // Secant
                case 4: // RegulaFalsi
                {
                    ImGui::Combo("Test function f(x)", &testFuncIdx, kTestFuncNames, IM_ARRAYSIZE(kTestFuncNames));
                    if (currentFunction == 2)
                        ImGui::TextDisabled("(Starts iterating from x = 1)");
                    else
                        ImGui::TextDisabled("(Automatically searches [-100, 100] for a sign change)");
                    if (ImGui::Button("Calculate"))
                    {
                        RootForge::Function f = GetTestFunction(testFuncIdx);
                        double r = 0.0; string name;
                        if (currentFunction == 2) { r = RootForge::Root::NewtonRaphson(f); name = "NewtonRaphson"; }
                        else if (currentFunction == 3) { r = RootForge::Root::Secant(f); name = "Secant"; }
                        else { r = RootForge::Root::RegulaFalsi(f); name = "RegulaFalsi"; }
                        ostringstream oss; oss << name << " root of " << kTestFuncNames[testFuncIdx] << " = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
            }
        }
        // ------------------------------------------------------------------
        // Library 4: Interpolation
        // ------------------------------------------------------------------
        else if (currentLibrary == 4)
        {
            switch (currentFunction)
            {
                case 0: // LinearInterpolation
                {
                    ImGui::InputFloat("x1", &valA);
                    ImGui::InputFloat("y1", &valB);
                    ImGui::InputFloat("x2", &valC);
                    ImGui::InputFloat("y2", &valD);
                    ImGui::InputFloat("x_i (point to interpolate)", &valE);
                    if (ImGui::Button("Calculate"))
                    {
                        double r = RootForge::Interpolation::LinearInterpolation(valA, valB, valC, valD, valE);
                        ostringstream oss; oss << "LinearInterpolation at x=" << valE << " = " << r;
                        resultText = oss.str(); history.push_back(resultText);
                    }
                    break;
                }
                case 1: // NumericalPolyInterpolation
                {
                    ImGui::InputText("x values (comma separated)", vecBufX, IM_ARRAYSIZE(vecBufX));
                    ImGui::InputText("y values (comma separated)", vecBufY, IM_ARRAYSIZE(vecBufY));
                    ImGui::InputFloat("x_i (point to interpolate)", &valE);
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> xs = ParseDoubles(vecBufX);
                        vector<double> ys = ParseDoubles(vecBufY);
                        if (xs.size() != ys.size() || xs.empty())
                        {
                            resultText = "x and y lists must be the same non-zero length.";
                        }
                        else
                        {
                            double r = RootForge::Interpolation::NumericalPolyInterpolation(xs, ys, valE);
                            ostringstream oss; oss << "Interpolated value at x=" << valE << " = " << r;
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 2: // LagrangePolyInterpolation
                {
                    ImGui::InputText("x values (comma separated)", vecBufX, IM_ARRAYSIZE(vecBufX));
                    ImGui::InputText("y values (comma separated)", vecBufY, IM_ARRAYSIZE(vecBufY));
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> xs = ParseDoubles(vecBufX);
                        vector<double> ys = ParseDoubles(vecBufY);
                        if (xs.size() != ys.size() || xs.empty())
                        {
                            resultText = "x and y lists must be the same non-zero length.";
                        }
                        else
                        {
                            vector<double> coeffs = RootForge::Interpolation::LagrangePolyInterpolation(xs, ys);
                            resultText = "Lagrange polynomial coefficients (low to high degree):\n" + FormatVec(coeffs);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
            }
        }
        // ------------------------------------------------------------------
        // Library 5: ODE Solvers
        // ------------------------------------------------------------------
        else if (currentLibrary == 5)
        {
            if (currentFunction != 4)
            {
                ImGui::Combo("dy/dx = f(x, y)", &odeFuncIdx, kODEFuncNames, IM_ARRAYSIZE(kODEFuncNames));
                ImGui::InputFloat("x0", &valA);
                ImGui::InputFloat("y0", &valB);
                ImGui::InputFloat("h (step size)", &valC);
                ImGui::InputFloat("target x", &valD);
                if (ImGui::Button("Calculate"))
                {
                    RootForge::DFunction f = GetODEFunction(odeFuncIdx);
                    double r = 0.0; string name;
                    switch (currentFunction)
                    {
                        case 0: r = RootForge::ODE::EulerMethod(f, valA, valB, valC, valD); name = "EulerMethod"; break;
                        case 1: r = RootForge::ODE::RK4(f, valA, valB, valC, valD); name = "RK4"; break;
                        case 2: r = RootForge::ODE::Midpoint(f, valA, valB, valC, valD); name = "Midpoint"; break;
                        case 3: r = RootForge::ODE::RK2(f, valA, valB, valC, valD); name = "RK2"; break;
                    }
                    ostringstream oss;
                    oss << name << ": y(" << valD << ") ~= " << r
                        << "   [" << kODEFuncNames[odeFuncIdx] << ", starting at (" << valA << ", " << valB << ")]";
                    resultText = oss.str(); history.push_back(resultText);
                }
            }
            else // SystemRK4 demo: y0' = y1, y1' = -y0 (simple harmonic oscillator)
            {
                ImGui::TextDisabled("Demo system: y0' = y1,  y1' = -y0  (simple harmonic oscillator)");
                ImGui::InputFloat("x0", &valA);
                ImGui::InputFloat("y0[0] (position)", &valB);
                ImGui::InputFloat("y0[1] (velocity)", &valC);
                ImGui::InputFloat("h (step size)", &valD);
                ImGui::InputFloat("target x", &valE);
                if (ImGui::Button("Calculate"))
                {
                    vector<RootForge::NFunction> system = {
                        [](double x, vector<double> y){ (void)x; return y[1]; },
                        [](double x, vector<double> y){ (void)x; return -y[0]; }
                    };
                    vector<double> y0 = { (double)valB, (double)valC };
                    vector<double> r = RootForge::ODE::SystemRK4(system, valA, y0, valD, valE);
                    ostringstream oss;
                    oss << "SystemRK4 at x=" << valE << ": y0 = " << r[0] << ", y1 = " << r[1];
                    resultText = oss.str(); history.push_back(resultText);
                }
            }
        }
        // ------------------------------------------------------------------
        // Library 6: Polynomial
        // ------------------------------------------------------------------
        else if (currentLibrary == 6)
        {
            switch (currentFunction)
            {
                case 0: // Evaluate
                {
                    ImGui::InputText("coefficients (low -> high degree)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter at least one coefficient."; }
                        else
                        {
                            Polynomial p(c);
                            double r = p.evaluate(valA);
                            ostringstream oss; oss << "p(" << valA << ") = " << r;
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 1: // Add
                case 2: // Subtract
                case 4: // Divide
                case 5: // Multiply
                {
                    ImGui::InputText("p1 coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::InputText("p2 coefficients (low -> high)", polyBufB, IM_ARRAYSIZE(polyBufB));
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c1 = ParseDoubles(polyBufA);
                        vector<double> c2 = ParseDoubles(polyBufB);
                        if (c1.empty() || c2.empty())
                        {
                            resultText = "Enter coefficients for both polynomials.";
                        }
                        else
                        {
                            Polynomial p1(c1), p2(c2);
                            if (currentFunction == 1)
                            {
                                Polynomial r = Polynomial::add(p1, p2);
                                resultText = "p1 + p2 coefficients:\n" + FormatVec(r.coeffs);
                            }
                            else if (currentFunction == 2)
                            {
                                Polynomial r = Polynomial::negate(p1, p2);
                                resultText = "p1 - p2 coefficients:\n" + FormatVec(r.coeffs);
                            }
                            else if (currentFunction == 4)
                            {
                                vector<Polynomial> r = Polynomial::divide(p1, p2);
                                resultText = "Quotient coefficients:\n" + FormatVec(r[0].coeffs) +
                                              "\nRemainder coefficients:\n" + FormatVec(r[1].coeffs);
                            }
                            else // multiply
                            {
                                vector<double> r = Polynomial::multiply(p1, p2);
                                resultText = "p1 * p2 coefficients:\n" + FormatVec(r);
                            }
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 3: // ScalarProd
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::InputFloat("scalar k", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter coefficients."; }
                        else
                        {
                            Polynomial p(c);
                            Polynomial r = Polynomial::ScalarProd(p, valA);
                            resultText = "k * p coefficients:\n" + FormatVec(r.coeffs);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 6: // NewtonHorner
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::TextDisabled("(Works best when all roots of p are real)");
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.size() < 2) { resultText = "Enter at least 2 coefficients (degree >= 1)."; }
                        else
                        {
                            Polynomial p(c);
                            vector<double> roots = Polynomial::NewtonHorner(p);
                            resultText = "Roots found:\n" + FormatVec(roots);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 7: // Differentiate
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter coefficients."; }
                        else
                        {
                            Polynomial p(c);
                            Polynomial r = Polynomial::Differentiate(p);
                            resultText = "p'(x) coefficients:\n" + FormatVec(r.coeffs);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 8: // Integrate
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter coefficients."; }
                        else
                        {
                            Polynomial p(c);
                            Polynomial r = Polynomial::Integrate(p);
                            resultText = "Antiderivative coefficients (C=0):\n" + FormatVec(r.coeffs);
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 9: // NumericalDiff
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::InputFloat("x", &valA);
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter coefficients."; }
                        else
                        {
                            Polynomial p(c);
                            double r = p.NumericalDiff(valA);
                            ostringstream oss; oss << "p'(" << valA << ") ~= " << r;
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 10: // NumericalInt
                {
                    ImGui::InputText("p coefficients (low -> high)", polyBufA, IM_ARRAYSIZE(polyBufA));
                    ImGui::InputFloat("x start", &valA);
                    ImGui::InputFloat("x end", &valB);
                    if (ImGui::Button("Calculate"))
                    {
                        vector<double> c = ParseDoubles(polyBufA);
                        if (c.empty()) { resultText = "Enter coefficients."; }
                        else
                        {
                            Polynomial p(c);
                            double r = p.NumericalInt(valA, valB);
                            ostringstream oss; oss << "Integral of p from " << valA << " to " << valB << " ~= " << r;
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
                case 11: // Quadratic Formula
                {
                    ImGui::InputFloat("a", &valA);
                    ImGui::InputFloat("b", &valB);
                    ImGui::InputFloat("c", &valC);
                    if (ImGui::Button("Calculate"))
                    {
                        double a = valA, b = valB, c = valC;
                        if (a == 0.0)
                        {
                            resultText = "'a' cannot be 0.";
                        }
                        else
                        {
                            double disc = b * b - 4 * a * c;
                            ostringstream oss; oss << fixed << setprecision(6);
                            if (disc > 0)
                            {
                                double r1 = (-b - sqrt(disc)) / (2 * a);
                                double r2 = (-b + sqrt(disc)) / (2 * a);
                                oss << "Discriminant = " << disc << " (two real roots)\nRoot 1 = " << r1 << "\nRoot 2 = " << r2;
                            }
                            else if (disc == 0)
                            {
                                double r = -b / (2 * a);
                                oss << "Discriminant = 0 (one repeated real root)\nRoot = " << r;
                            }
                            else
                            {
                                double realPart = -b / (2 * a);
                                double imagVal = (-disc) / (4 * a * a);
                                string imagStr = RootForge::MathematicalFunctions::simplify_sqrt((int)round(imagVal));
                                oss << "Discriminant = " << disc << " (complex roots)\nRoot 1 = " << realPart << " + " << imagStr
                                    << "i\nRoot 2 = " << realPart << " - " << imagStr << "i";
                            }
                            // Also print the library's own console report:
                            Polynomial::quadratic_formula(a, b, c);
                            resultText = oss.str();
                        }
                        history.push_back(resultText);
                    }
                    break;
                }
            }
        }

        // --- Result panel ---
        ImGui::Separator();
        ImGui::Text("Result:");
        ImGui::TextWrapped("%s", resultText.c_str());

        // --- History panel ---
        ImGui::Separator();
        if (ImGui::Button("Clear History")) history.clear();
        ImGui::BeginChild("HistoryRegion", ImVec2(0, 180), true);
        for (int i = (int)history.size() - 1; i >= 0; i--)
        {
            ImGui::TextWrapped("%s", history[i].c_str());
            ImGui::Separator();
        }
        ImGui::EndChild();

        ImGui::End();
        // --- YOUR UI LOGIC ENDS HERE ---

        // 5. Render to the Screen
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // 6. Cleanup and Shutdown
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
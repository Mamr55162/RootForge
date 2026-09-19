#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <cstddef>
#include <vector>
using namespace std;

#ifndef ROOTFORGE_ROOTFORGE_H
#define ROOTFORGE_ROOTFORGE_H

#endif //ROOTFORGE_ROOTFORGE_H


namespace RootForge
{
    using Function = std::function<double(double)>;
    using DFunction = std::function<double(double, double)>;
    using NFunction = std::function<double(double, vector<double>)>;
    enum class Status
    {
        Success,
        InvalidInput,
        DivisionByZero,
        MaxIterations,
        NotConverged
    };

    struct NumericalResult
    {
        double value;
        double absError;
        double relativeError;
        double percentError;
        std::size_t iterations;
        Status status;

        static void ShowResult(const NumericalResult result)
        {
            cout << "Result = " << result.value << endl;
            cout << "Absolute Error = " << result.absError << endl;
            cout << "Relative Error = " << result.relativeError << endl;
            cout << "Percentage Error = " << result.percentError << " %" << endl;
            cout << "Iterations = " << result.iterations << endl;
        }
    };
}

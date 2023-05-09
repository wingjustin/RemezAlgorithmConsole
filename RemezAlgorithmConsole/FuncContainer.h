#pragma once

#include "DoubleFunc.h"

namespace MyMath {
    struct FuncContainer : DoubleFunc {
        double (*func)(double);

        FuncContainer(double (*func)(double));

        double operator()(double x);
    };

    struct EstimateFuncContainer : DoubleFunc {
        unsigned int numeratorCount;
        unsigned int denominatorCount;

        bool pinned;

        double* coefficients;

        EstimateFuncContainer(double* coefficients, bool pinnied
            , unsigned int numeratorCount, unsigned int denominatorCount);

        double operator()(double x);
	};
}

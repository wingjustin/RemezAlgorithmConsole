#pragma once

#include "DoubleFunc.h"

namespace MyMath {
    struct FuncContainer : DoubleFunc {
        double (*func)(double);

        FuncContainer(double (*func)(double));

        double operator()(double x);
    };

    struct RationalFuncContainer : DoubleFunc {
        unsigned int numeratorCount;
        unsigned int denominatorCount;

        bool pinned;

        double* coefficients;

        RationalFuncContainer(double* coefficients, bool pinnied
            , unsigned int numeratorCount, unsigned int denominatorCount);

        double operator()(double x);
	};
}

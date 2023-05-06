#pragma once

#include "FuncContainer.h"

#ifndef _DOUBLE_HI
#define _DOUBLE_HI(x) *(1 + (int*)&x)
#endif

namespace MyMath {
    struct RelativeErrorFunc : DoubleFunc {
        RationalFuncContainer* rationalFunc;
        FuncContainer* oriFunc;

        RelativeErrorFunc(RationalFuncContainer* rationalFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct AbsoluteErrorFunc : DoubleFunc {
        RationalFuncContainer* rationalFunc;
        FuncContainer* oriFunc;

        AbsoluteErrorFunc(RationalFuncContainer* rationalFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct Abs_RelativeErrorFunc : RelativeErrorFunc {

        Abs_RelativeErrorFunc(RationalFuncContainer* rationalFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct Abs_AbsoluteErrorFunc : AbsoluteErrorFunc {

        Abs_AbsoluteErrorFunc(RationalFuncContainer* rationalFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };
}

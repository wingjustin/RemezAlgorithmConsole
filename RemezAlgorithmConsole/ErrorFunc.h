#pragma once

#include "FuncContainer.h"

#ifndef _DOUBLE_HI
#define _DOUBLE_HI(x) *(1 + (int*)&x)
#endif

namespace MyMath {
    struct WeightedErrorFunc : DoubleFunc {
        EstimateFuncContainer* estimFunc;
        FuncContainer* oriFunc;
        FuncContainer* weightFunc;

        WeightedErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc, FuncContainer* weightFunc);

        double operator()(double x);
    };

    struct RelativeErrorFunc : DoubleFunc {
        EstimateFuncContainer* estimFunc;
        FuncContainer* oriFunc;

        RelativeErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct AbsoluteErrorFunc : DoubleFunc {
        EstimateFuncContainer* estimFunc;
        FuncContainer* oriFunc;

        AbsoluteErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct Abs_WeightedErrorFunc : WeightedErrorFunc {

        Abs_WeightedErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc, FuncContainer* weightFunc);

        double operator()(double x);
    };

    struct Abs_RelativeErrorFunc : RelativeErrorFunc {

        Abs_RelativeErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };

    struct Abs_AbsoluteErrorFunc : AbsoluteErrorFunc {

        Abs_AbsoluteErrorFunc(EstimateFuncContainer* estimFunc, FuncContainer* oriFunc);

        double operator()(double x);
    };
}

#include "ErrorFunc.h"

using namespace MyMath;

WeightedErrorFunc::WeightedErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc, FuncContainer* const weightFunc) {
    WeightedErrorFunc::estimFunc = estimFunc;
    WeightedErrorFunc::oriFunc = oriFunc;
    WeightedErrorFunc::weightFunc = weightFunc;
}

double WeightedErrorFunc::operator()(const double x) {
    double w = (*weightFunc)(x);
    return w
        ? ((*oriFunc)(x) - (*estimFunc)(x)) / w // (f-p)/w
        : (*oriFunc)(x) - (*estimFunc)(x); // (f-p)/w
}

RelativeErrorFunc::RelativeErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc) {
    RelativeErrorFunc::estimFunc = estimFunc;
    RelativeErrorFunc::oriFunc = oriFunc;
}

double RelativeErrorFunc::operator()(const double x) {
    double f // true value
        , err = (f = (*oriFunc)(x)) - (*estimFunc)(x); // f-p
    _DOUBLE_HI(f) &= 0x7FFFFFFF; // get absolute, |f|

    return f ? err / f : err;
}

AbsoluteErrorFunc::AbsoluteErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc) {
    AbsoluteErrorFunc::estimFunc = estimFunc;
    AbsoluteErrorFunc::oriFunc = oriFunc;
}

double AbsoluteErrorFunc::operator()(const double x) {
    return (*oriFunc)(x) - (*estimFunc)(x); // f-p
}

Abs_WeightedErrorFunc::Abs_WeightedErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc, FuncContainer* const weightFunc)
    :WeightedErrorFunc(estimFunc, oriFunc, weightFunc) {
}

double Abs_WeightedErrorFunc::operator()(const double x) {
    double absErr = WeightedErrorFunc::operator()(x); // (f-p)/w
    _DOUBLE_HI(absErr) &= 0x7FFFFFFF;
    return absErr;
}

Abs_RelativeErrorFunc::Abs_RelativeErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc)
    :RelativeErrorFunc(estimFunc, oriFunc) {
}

double Abs_RelativeErrorFunc::operator()(const double x) {
    double absErr = RelativeErrorFunc::operator()(x); // (f-p)/|f|
    _DOUBLE_HI(absErr) &= 0x7FFFFFFF;
    return absErr;
}

Abs_AbsoluteErrorFunc::Abs_AbsoluteErrorFunc(EstimateFuncContainer* const estimFunc, FuncContainer* const oriFunc)
    :AbsoluteErrorFunc(estimFunc, oriFunc) {
}

double Abs_AbsoluteErrorFunc::operator()(const double x) {
    double absErr = AbsoluteErrorFunc::operator()(x); // f-p
    _DOUBLE_HI(absErr) &= 0x7FFFFFFF;
    return absErr;
}

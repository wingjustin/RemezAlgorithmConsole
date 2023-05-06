#include "ErrorFunc.h"

using namespace MyMath;

RelativeErrorFunc::RelativeErrorFunc(RationalFuncContainer* const rationalFunc, FuncContainer* const oriFunc) {
    RelativeErrorFunc::rationalFunc = rationalFunc;
    RelativeErrorFunc::oriFunc = oriFunc;
}

double RelativeErrorFunc::operator()(const double x) {
    double f // true value
        , err = (f = (*oriFunc)(x)) - (*rationalFunc)(x); // f-p
    _DOUBLE_HI(f) &= 0x7FFFFFFF; // get absolute, |f|

    return f ? err / f : err;
}

AbsoluteErrorFunc::AbsoluteErrorFunc(RationalFuncContainer* const rationalFunc, FuncContainer* const oriFunc) {
    AbsoluteErrorFunc::rationalFunc = rationalFunc;
    AbsoluteErrorFunc::oriFunc = oriFunc;
}

double AbsoluteErrorFunc::operator()(const double x) {
    return (*oriFunc)(x) - (*rationalFunc)(x); // f-p
}

Abs_RelativeErrorFunc::Abs_RelativeErrorFunc(RationalFuncContainer* const rationalFunc, FuncContainer* const oriFunc)
    :RelativeErrorFunc(rationalFunc, oriFunc) {
}

double Abs_RelativeErrorFunc::operator()(const double x) {
    double absErr = RelativeErrorFunc::operator()(x); // (f-p)/|f|
    _DOUBLE_HI(absErr) &= 0x7FFFFFFF;
    return absErr;
}

Abs_AbsoluteErrorFunc::Abs_AbsoluteErrorFunc(RationalFuncContainer* const rationalFunc, FuncContainer* const oriFunc)
    :AbsoluteErrorFunc(rationalFunc, oriFunc) {
}

double Abs_AbsoluteErrorFunc::operator()(const double x) {
    double absErr = AbsoluteErrorFunc::operator()(x); // f-p
    _DOUBLE_HI(absErr) &= 0x7FFFFFFF;
    return absErr;
}

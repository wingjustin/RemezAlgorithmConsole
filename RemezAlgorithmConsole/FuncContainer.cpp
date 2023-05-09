#include "FuncContainer.h"

using namespace MyMath;

EstimateFuncContainer::EstimateFuncContainer(double* const coefficients, const bool pinnied
    , const unsigned int numeratorCount, const unsigned int denominatorCount) {
    EstimateFuncContainer::coefficients = coefficients;
    EstimateFuncContainer::pinned = pinnied;
    EstimateFuncContainer::numeratorCount = numeratorCount;
    EstimateFuncContainer::denominatorCount = denominatorCount;
}

double EstimateFuncContainer::operator()(const double x) {
    double n = 0.0; // numerator
    double d = 0.0; // denominator

    //calculate numerator
    switch (numeratorCount) {
    default:
        n = coefficients[numeratorCount - 1];
        for (unsigned int i = numeratorCount - 2; i > 0; i--) {
            n *= x;
            n += coefficients[i];
        }
        n *= x;
    case 1:
        n += coefficients[0];
        if (pinned)
            n *= x;
    case 0:
        break;
    }

    //calculate denominator
    switch (denominatorCount) {
    default:
        d = coefficients[numeratorCount + denominatorCount - 1];
        for (unsigned int i = numeratorCount + denominatorCount - 2; i > numeratorCount; i--) {
            d *= x;
            d += coefficients[i];
        }
        d *= x;
    case 1:
        d += coefficients[numeratorCount];
        d *= x;
    case 0:
        d++;
        break;
    }

    return n / d; // approximate
}

FuncContainer::FuncContainer(double (* const func)(double)) {
    FuncContainer::func = func;
}

double FuncContainer::operator()(const double x) {
    return func(x);
}

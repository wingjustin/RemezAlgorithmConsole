#pragma once

#include "DoubleFunc.h"

#ifndef _DOUBLE_HI
#define _DOUBLE_HI(x) *(1 + (int*)&x)
#endif

#ifndef _BRENTMETHOD_GOLDEN_RATIO
#define _BRENTMETHOD_GOLDEN_RATIO 1.6180339887498948482045868343656
#define _BRENTMETHOD_GOLDEN_RATIO_A 0.61803398874989484820458683436564
#define _BRENTMETHOD_GOLDEN_RATIO_B 0.38196601125010515179541316563436
#endif

#ifndef _BRENTMETHOD_SQRT_DOUBLE_EPSILON
#define _BRENTMETHOD_SQRT_DOUBLE_EPSILON 0.00000001490116119384765625 //std::sqrt(std::numeric_limits<double>::epsilon())
#define _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER 0.0000000037252902984619140625 //quarter of _SQRT_DOUBLE_EPSILON
#endif

namespace MyMath {
	class BrentMethodExtrema{
	public:
		static double* FindMinima(double (*func)(double), double left, double right, unsigned int& max_iter); // return x
		static double* FindMaxima(double (*func)(double), double left, double right, unsigned int& max_iter); // return x

		static double* FindMinima(DoubleFunc* func, double left, double right, unsigned int& max_iter); // return x
		static double* FindMaxima(DoubleFunc* func, double left, double right, unsigned int& max_iter); // return x
	};
}

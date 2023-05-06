#pragma once

#include "DoubleFunc.h"

#ifndef _DOUBLE_HI
#define _DOUBLE_HI(x) *(1 + (int*)&x)
#endif

#ifndef _TOMS748_DOUBLE_EPSILON
#define _TOMS748_DOUBLE_EPSILON 2.2204460492503131e-016 //std::numeric_limits<double>::epsilon()
#define _TOMS748_DOUBLE_DOUBLE_EPSILON 4.4408920985006262e-016 // 2.0 * _TOMS748_DOUBLE_EPSILON
#define _TOMS748_FOUR_TIMES_DOUBLE_EPSILON 8.8817841970012524e-016 // 4.0 * _TOMS748_DOUBLE_EPSILON
#define _TOMS748_FIVE_TIMES_DOUBLE_EPSILON 1.11022302462515655e-015 // 5.0 * _TOMS748_DOUBLE_EPSILON

#define _TOMS748_DOUBLE_MIN 2.2250738585072014e-308 //std::numeric_limits<double>::min();
#define _TOMS748_32_TIMES_DOUBLE_MIN 7.12023634722304448e-307 // 32.0 * _TOMS748_DOUBLE_MIN

#define _TOMS748_MU 0.5
#endif

/****
* references:
* 1. boost.org, \boost\math\tools\toms748_solve.hpp
* 2. document : Alefeld, Potra and Shi: 1995 Algorithm 748 Enclosing Zeros of Continuous Functions
* 3. https://github.com/jamadagni/toms748/blob/master/toms748.cpp
*
* The signs of f(a) and f(b) must be different.
*
* There are four main methods
* 1. Bisection
* 2. Secant interpolation, and Double-Length Secant Interpolate
* 3. Newton Quadratic interpolation, (just like Newton Polynomial + Newton's Method)
* 4. Inverse Cubic interpolation, (just like Cubic version of Inverse Quadratic interpolation in Brent-Dekker Method)
*/

namespace MyMath {
	class TOMS748Algorithm {
	private:
		static bool CheckValidNum(double& x);
		static bool PassTolerance(double& a, double& b);
		static bool MeetICIConditions(double& fa, double& fb, double& fd, double& fe);

		//methods
		static void UpdateBracket(double(&func)(double), double& c, double& a, double& b, double& d
			, double& fc, double& fa, double& fb, double& fd);
		static void UpdateBracket(DoubleFunc& func, double& c, double& a, double& b, double& d
			, double& fc, double& fa, double& fb, double& fd);

		static void SecantInterpolate(double& c, double& a, double& b
			, double& fa, double& fb);
		static void DoubleLengthSecantInterpolate(double& c, double& a, double& b
			, double& fa, double& fb);
		static void NewtonQuadraticInterpolate(unsigned int step, double& c, double& a, double& b, double& d
			, double& fa, double& fb, double& fd);
		static void InverseCubicInterpolate(double& c, double& a, double& b, double& d, double& e
			, double& fa, double& fb, double& fd, double& fe);
	public:
		static double FindRoot(double (*func)(double), double left, double right, unsigned int& max_iter); // return x
		static double FindRoot(DoubleFunc* func, double left, double right, unsigned int& max_iter); // return x
	};
}

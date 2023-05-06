#include "TOMS748Algorithm.h"

using namespace MyMath;

/****
* references:
* 1. boost.org, .\boost\math\tools\toms748_solve.hpp
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

bool TOMS748Algorithm::CheckValidNum(double& x) {
	return (_DOUBLE_HI(x) & 0x7FF00000) != 0x7FF00000; // check if exp bits are all 1, then it may be +inf, -inf or +-NaN.
}

bool TOMS748Algorithm::PassTolerance(double& a, double& b) {
	double abs_a = a, abs_b = b, abs_temp = a - b;
	_DOUBLE_HI(abs_a) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_b) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_temp) &= 0x7FFFFFFF;
	return abs_temp > (_TOMS748_FOUR_TIMES_DOUBLE_EPSILON * (abs_a < abs_b ? abs_a : abs_b)); // 4.0 * EPSILON * (|a| or |b|)
}

bool TOMS748Algorithm::MeetICIConditions(double& fa, double& fb, double& fd, double& fe) {
	//make sure all four function values fa, fb, fd, and fe are distinct

	//get absolute values
	double abs;

	abs = fa - fb;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	abs = fa - fd;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	abs = fa - fe;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	abs = fb - fd;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	abs = fb - fe;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	abs = fd - fe;
	_DOUBLE_HI(abs) &= 0x7FFFFFFF;
	if (abs < _TOMS748_32_TIMES_DOUBLE_MIN)
		return false;

	return true;
}

void TOMS748Algorithm::UpdateBracket(double(&func)(double), double& c, double& a, double& b, double& d
	, double& fc, double& fa, double& fb, double& fd) {
	/*
	* a is the second best approximation root
	* c is current guess
	* d is the third best guess, the last removed point, a or b
	*/

	//update bracket
	// check if c is too close to a or b
	double abs_a, abs_b, limit_a, limit_b;
	abs_a = a, abs_b = b;
	_DOUBLE_HI(abs_a) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_b) &= 0x7FFFFFFF;
	limit_a = abs_a * _TOMS748_DOUBLE_DOUBLE_EPSILON; // 2.0 * EPSILON
	limit_b = abs_b * _TOMS748_DOUBLE_DOUBLE_EPSILON; // 2.0 * EPSILON

	if ((b - a) < _TOMS748_FOUR_TIMES_DOUBLE_EPSILON * a) // 4.0 * EPSILON
		c = a + (b - a) * 0.5; //avoid large num + large num, (a + b) * 0.5; //Bisection
	else if (c <= a + limit_a)
		c = a + limit_a;
	else if (c >= b - limit_b)
		c = b - limit_b;

	fc = func(c);

	if (fc == 0) {
		a = c;
		fa = 0;
		d = 0;
		fd = 0;
		return;
	}

	//Bisection
	// check sign
	if ((_DOUBLE_HI(fc) & 0x80000000) == (_DOUBLE_HI(fa) & 0x80000000)) {
		//change the interval to [c, b]
		d = a;
		fd = fa;
		a = c;
		fa = fc;
	}
	else {
		//change the interval to [a, c]
		d = b;
		fd = fb;
		b = c;
		fb = fc;
	}
}

void TOMS748Algorithm::UpdateBracket(DoubleFunc& func, double& c, double& a, double& b, double& d
	, double& fc, double& fa, double& fb, double& fd) {
	/*
	* a is the second best approximation root
	* c is current guess
	* d is the third best guess, the last removed point, a or b
	*/

	//update bracket
	// check if c is too close to a or b
	double abs_a, abs_b, limit_a, limit_b;
	abs_a = a, abs_b = b;
	_DOUBLE_HI(abs_a) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_b) &= 0x7FFFFFFF;
	limit_a = abs_a * _TOMS748_DOUBLE_DOUBLE_EPSILON; // 2.0 * EPSILON
	limit_b = abs_b * _TOMS748_DOUBLE_DOUBLE_EPSILON; // 2.0 * EPSILON

	if ((b - a) < _TOMS748_FOUR_TIMES_DOUBLE_EPSILON * a) // 4.0 * EPSILON
		c = a + (b - a) * 0.5; //avoid large num + large num, (a + b) * 0.5; //Bisection
	else if (c <= a + limit_a)
		c = a + limit_a;
	else if (c >= b - limit_b)
		c = b - limit_b;

	fc = func(c);

	if (fc == 0) {
		a = c;
		fa = 0;
		d = 0;
		fd = 0;
		return;
	}

	//Bisection
	// check sign
	if ((_DOUBLE_HI(fc) & 0x80000000) == (_DOUBLE_HI(fa) & 0x80000000)) {
		//change the interval to [c, b]
		d = a;
		fd = fa;
		a = c;
		fa = fc;
	}
	else {
		//change the interval to [a, c]
		d = b;
		fd = fb;
		b = c;
		fb = fc;
	}
}

void TOMS748Algorithm::SecantInterpolate(double& c, double& a, double& b
	, double& fa, double& fb) {
	//Secant interpolation

	double abs_a, abs_b, limit_a, limit_b;
	abs_a = a, abs_b = b;
	_DOUBLE_HI(abs_a) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_b) &= 0x7FFFFFFF;
	limit_a = abs_a * _TOMS748_FIVE_TIMES_DOUBLE_EPSILON; // 5.0 * EPSILON
	limit_b = abs_b * _TOMS748_FIVE_TIMES_DOUBLE_EPSILON; // 5.0 * EPSILON

	c = fa / (fa - fb) * (b - a) + a;
	// check if c overflow
	if (c <= a + limit_a || c >= b - limit_b)
		c = a + (b - a) * 0.5; //avoid large num + large num, (a + b) * 0.5; //Bisection
}

void TOMS748Algorithm::DoubleLengthSecantInterpolate(double& c, double& a, double& b
	, double& fa, double& fb) {
	//Double-Length Secant interpolation

	double u, fu;
	double abs_fa, abs_fb, abs_cmu;
	abs_fa = fa, abs_fb = fb;
	_DOUBLE_HI(abs_fa) &= 0x7FFFFFFF;
	_DOUBLE_HI(abs_fb) &= 0x7FFFFFFF;

	if (abs_fa < abs_fb) {
		u = a;
		fu = fa;
	}
	else {
		u = b;
		fu = fb;
	}
	c = (fu / (fb - fa)) * (b - a);
	c = u - c - c;

	abs_cmu = c - u;
	_DOUBLE_HI(abs_cmu) &= 0x7FFFFFFF;
	u = (b - a) * 0.5;

	if (abs_cmu > u)
		c = a + u; // c = (a + b) * 0.5; Bisection
}

void TOMS748Algorithm::NewtonQuadraticInterpolate(const unsigned int step, double& c, double& a, double& b, double& d
	, double& fa, double& fb, double& fd) {
	//Newton Quadratic interpolation
	// Newton Quadratic Polynomial / Newton's divided differences interpolation polynomial
	// P(x) = P(a, b, d)(x)
	// P(x) = f(a) + f[a, b](x - a) + f[a, b, d](x - a)(x - b)
	// Newton's Method
	// x1 = x0 - f(x0)/f'(x0)
	// x2 = x1 - f(x1)/f'(x1)
	// x3 = x2 - f(x2)/f'(x2)
	// ...
	// let A = f[a, b, d], B = f[a, b]
	// P(x)
	// = f(a) + B(x - a) + A(x - a)(x - b)
	// = f(a) + B*x - B*a + A(x*x - ax - bx + a*b)
	// = A*x*x + (B - A*a - A*b)*x + (A*a*b - B*a + f(a))
	// so, P'(x) = 2*A*x + (B - A*a - A*b)
	//  = A*(2*x - a - b) + B
	// P(x) = f(a) + B(x - a) + A(x - a)(x - b)
	//  = f(a) + (x - a)*[B + A(x - b)]
	// next_c = c - P(c) / P'(c)

	double A, B;

	B = (fb - fa) / (b - a); // B = f[a, b]
	if(!CheckValidNum(B))
		return SecantInterpolate(c, a, b, fa, fb);

	A = (fd - fb) / (d - b); // A = f[a, b, d] = (f[b, d] - f[a, b]) / (d - a)
	if (!CheckValidNum(A))
		return SecantInterpolate(c, a, b, fa, fb);

	A = (A - B) / (d - a);
	if (A == 0 || !CheckValidNum(A))
		return SecantInterpolate(c, a, b, fa, fb);

	// Determine the starting point of the Newton steps:
	// A is coefficient of x**2 in P(x)
	// so the parabola opens up when A is positive, opens down when A is negative
	// check sign
	c = (_DOUBLE_HI(A) & 0x80000000) == (_DOUBLE_HI(fa) & 0x80000000)
		? a
		: b;

	// Take the Newton steps:
	double next_c;
	unsigned int i = 0;
	for (; i < step; i++) {
		next_c = c
			- (fa + (c - a) * (B + A * (c - b)))
			/ (A * (2 * c - a - b) + B);

		if (!CheckValidNum(next_c))
			next_c = a - 1; // falls out the bracket

		// if it falls out of bracket, get the result of previous step 
		if ((next_c <= a) || (next_c >= b))
			break;
		else
			c = next_c;
	}
	//check if any Newton steps are success
	if (i == 0)
		return SecantInterpolate(c, a, b, fa, fb); // no one step is success
}

void TOMS748Algorithm::InverseCubicInterpolate(double& c, double& a, double& b, double& d, double& e
	, double& fa, double& fb, double& fd, double& fe) {
	// formula references:
	// 1. boost.org, .\boost\math\tools\toms748_solve.hpp
	// 2. document : Alefeld, Potra and Shi: 1995 Algorithm 748 Enclosing Zeros of Continuous Functions, pages 332 to 334
	//
	// Points d and e lie outside the interval [a, b]
	// d is the third best guess
	// e is the forth best guess, the last value of d

	double q11 = (d - e) * fd / (fe - fd);
	double q21 = (b - d) * fb / (fd - fb);
	double q31 = (a - b) * fa / (fb - fa);
	double d21 = (b - d) * fd / (fd - fb);
	double d31 = (a - b) * fb / (fb - fa);

	double q22 = (d21 - q11) * fb / (fe - fb);
	double q32 = (d31 - q21) * fa / (fd - fa);
	double d32 = (d31 - q21) * fd / (fd - fa);
	double q33 = (d32 - q22) * fa / (fe - fa);
	c = q31 + q32 + q33 + a;

	// check if c is overflow
	if ((c <= a) || (c >= b))
		return NewtonQuadraticInterpolate(3, c, a, b, d, fa, fb, fd);
}

double TOMS748Algorithm::FindRoot(double (* const func)(double), const double left, const double right
	, unsigned int& max_iter) { // return x
	/*
	* a is the second best approximation root
	* c is current guess
	* d is the third best guess, the last removed point, a or b
	* e is the forth best guess, the last value of d
	* Points d and e lie outside the interval [a, b]
	*/

	if (max_iter == 0) // special case
		return left;

	unsigned int count = max_iter;

	double a, b, fa, fb;
	if (left < right) {
		a = left;
		b = right;
	}
	else if (left > right) {
		a = right;
		b = left;
	}
	else// if (left == right)
		return left; // cannot handle

	if (!PassTolerance(a, b)) {
		max_iter = 0;
		return a;
		//return b;
	}

	fa = func(a);
	fb = func(b);

	if (fa == 0) {
		max_iter = 0;
		return a;
	}
	if (fb == 0) {
		max_iter = 0;
		return b;
	}

	// check sign
	if ((_DOUBLE_HI(fa) & 0x80000000) == (_DOUBLE_HI(fb) & 0x80000000)) {
		max_iter = 0;
		return a;
		//return b;
	}

	double c, fc;
	double d, e, fd, fe;

	//fe = e = fd = 100000;

	// On the first step we take a secant step:
	SecantInterpolate(c, a, b, fa, fb);
	UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd); // gain d & fd
	--count;

	if (count && (fa != 0) && PassTolerance(a, b))
	{
		// On the second step we take a quadratic interpolation:
		NewtonQuadraticInterpolate(2, c, a, b, d, fa, fb, fd);
		e = d; //gain e & fe, the last value of d
		fe = fd;
		UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);
		--count;

		double initialDistance;

		while (count && (fa != 0) && PassTolerance(a, b))
		{
			// save initial distance
			initialDistance = b - a;

			//
			// Starting with the third step taken
			// we can use either quadratic or cubic interpolation.
			// Cubic interpolation requires that all four function values
			// fa, fb, fd, and fe are distinct, should that not be the case
			// then variable prof will get set to true, and we'll end up
			// taking a quadratic step instead.
			// boost.org
			//

			if (MeetICIConditions(fa, fb, fd, fe))
				InverseCubicInterpolate(c, a, b, d, e, fa, fb, fd, fe);
			else
				NewtonQuadraticInterpolate(2, c, a, b, d, fa, fb, fd);

			// bracket
			e = d;
			fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			//check for termination
			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			// another interpolated step:
			if (MeetICIConditions(fa, fb, fd, fe))
				InverseCubicInterpolate(c, a, b, d, e, fa, fb, fd, fe);
			else
				NewtonQuadraticInterpolate(3, c, a, b, d, fa, fb, fd);

			// Bracket again, and check termination condition
			//e = d;
			//fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			//double-length secant step:
			DoubleLengthSecantInterpolate(c, a, b, fa, fb);

			// Bracket again, and check termination condition:
			e = d;
			fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			//check for termination
			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			// check if the converging is not enough fast, the do an additional bisection step
			if ((b - a) >= _TOMS748_MU * initialDistance) {
				// bracket again
				e = d;
				fe = fd;
				c = a + (b - a) * 0.5; //avoid large num + large num, (a + b) * 0.5; //Bisection
				UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);
				--count;
			}
		} // while loop
	}

	max_iter -= count;
	if (fa == 0)
		return a;
	else if (fb == 0)
		return b;

	return a;
}

double TOMS748Algorithm::FindRoot(DoubleFunc* const func, const double left, const double right
	, unsigned int& max_iter) { // return x
	/*
	* a is the second best approximation root
	* c is current guess
	* d is the third best guess, the last removed point, a or b
	* e is the forth best guess, the last value of d
	* Points d and e lie outside the interval [a, b]
	*/

	if (max_iter == 0) // special case
		return left;

	unsigned int count = max_iter;

	double a, b, fa, fb;
	if (left < right) {
		a = left;
		b = right;
	}
	else if (left > right) {
		a = right;
		b = left;
	}
	else// if (left == right)
		return left; // cannot handle

	if (!PassTolerance(a, b)) {
		max_iter = 0;
		return a;
		//return b;
	}

	fa = (*func)(a);
	fb = (*func)(b);

	if (fa == 0) {
		max_iter = 0;
		return a;
	}
	if (fb == 0) {
		max_iter = 0;
		return b;
	}

	// check sign
	if ((_DOUBLE_HI(fa) & 0x80000000) == (_DOUBLE_HI(fb) & 0x80000000)) {
		max_iter = 0;
		return a;
		//return b;
	}

	double c, fc;
	double d, e, fd, fe;

	//fe = e = fd = 100000;

	// On the first step we take a secant step:
	SecantInterpolate(c, a, b, fa, fb);
	UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd); // gain d & fd
	--count;

	if (count && (fa != 0) && PassTolerance(a, b))
	{
		// On the second step we take a quadratic interpolation:
		NewtonQuadraticInterpolate(2, c, a, b, d, fa, fb, fd);
		e = d; //gain e & fe, the last value of d
		fe = fd;
		UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);
		--count;

		double initialDistance;

		while (count && (fa != 0) && PassTolerance(a, b))
		{
			// save initial distance
			initialDistance = b - a;

			//
			// Starting with the third step taken
			// we can use either quadratic or cubic interpolation.
			// Cubic interpolation requires that all four function values
			// fa, fb, fd, and fe are distinct, should that not be the case
			// then variable prof will get set to true, and we'll end up
			// taking a quadratic step instead.
			// boost.org
			//

			if (MeetICIConditions(fa, fb, fd, fe))
				InverseCubicInterpolate(c, a, b, d, e, fa, fb, fd, fe);
			else
				NewtonQuadraticInterpolate(2, c, a, b, d, fa, fb, fd);

			// bracket
			e = d;
			fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			//check for termination
			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			// another interpolated step:
			if (MeetICIConditions(fa, fb, fd, fe))
				InverseCubicInterpolate(c, a, b, d, e, fa, fb, fd, fe);
			else
				NewtonQuadraticInterpolate(3, c, a, b, d, fa, fb, fd);

			// Bracket again, and check termination condition
			//e = d;
			//fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			//double-length secant step:
			DoubleLengthSecantInterpolate(c, a, b, fa, fb);

			// Bracket again, and check termination condition:
			e = d;
			fe = fd;
			UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);

			//check for termination
			if (!(--count && (fa != 0) && PassTolerance(a, b)))
				break;

			// check if the converging is not enough fast, the do an additional bisection step
			if ((b - a) >= _TOMS748_MU * initialDistance) {
				// bracket again
				e = d;
				fe = fd;
				c = a + (b - a) * 0.5; //avoid large num + large num, (a + b) * 0.5; //Bisection
				UpdateBracket(*func, c, a, b, d, fc, fa, fb, fd);
				--count;
			}
		} // while loop
	}

	max_iter -= count;
	if (fa == 0)
		return a;
	else if (fb == 0)
		return b;

	return a;
}

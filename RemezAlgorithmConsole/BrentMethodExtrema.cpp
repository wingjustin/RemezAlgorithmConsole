#include "BrentMethodExtrema.h"

using namespace MyMath;

double* BrentMethodExtrema::FindMinima(double (* const func)(double), double const left, double const right, unsigned int& max_iter) {
	//unsigned int iters = 0;
	//unsigned int total_iters = 0;
	//unsigned int golden_sec_iters = 0;
	//unsigned int parabolic_iters = 0;

	/**
	* x is the least value of function
	* w is the second lowest value of function
	* v is previous value of w. usually it is next the lowest value
	*
	* u is last evaluated point
	* a is left bound
	* b is right bound
	* m is interval midpoint
	*
	* tol is tolerance
	* tol2 is tolerance2
	*
	* -----------------------------------------------------------------/
	** /--parabolic interpolation variables--/
	* p is numerator
	* q is denominator
	* u = x + p/q
	*
	* d is distance between x and u, also equal p/q. The distance moved in the last step
	* e is the second last p/q. The distance moved in the step before last
	*
	** /--conditions--/
	* q != 0;
	* u belong to [a,b]
	* abs(p/q) < 0.5 * abs(e)
	* pass tolerance check
	*
	* condition "p/q < n" can change to "p < n*q", when q > 0.
	* changing to "p < n*q" can also avoid q == 0
	*/

	if (max_iter == 0) // special case
		return new double[2] {left, func(left)}; //return left;

	unsigned int count = max_iter;

	double a, b, m;
	if (left < right) {
		a = left;
		b = right;
	}
	else {
		a = right;
		b = left;
	}

	double x, w, v, u, fx, fw, fv, fu;

	v = w = x = a + (b - a) * _BRENTMETHOD_GOLDEN_RATIO_B;
	fv = fw = fx = func(x);

	double d, e;
	d = e = 0;

	double r, p, q;
	double absP, absQ;
	double absTemp;

	double tol, tol2;
	// distance between x and u, check tol
	// distance between x and boundaries, check tol2

	//update paramters
	m = a + (b - a) * 0.5;
	absTemp = x;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
	tol2 = tol + tol;

	absTemp = x - m;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;

	// f must not be evaluated too close to a or b
	while (count && absTemp > (tol2 - (b - a) * 0.5)) {
		//iters++;
		//check SPI behaving
		// p = (w - x) * (w - x) * (fx - fv) + (v - x) * (v - x) * (fw - fx);
		// q = (w - x) * (fx - fv) + (v - x) * (fw - fx);
		// u = x + p / (q + q);
		absTemp = e;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp > tol) {
			// fit parabolic
			r = (v - x) * (fw - fx);
			q = (w - x) * (fx - fv);
			p = (w - x) * q;
			q += r;
			p += (v - x) * r;
			//q += q;
			r = e; // r is the second last p/2q
			e = d; // e becomes last p/2q

			// if(abs(p/(2*q)) < 0.5 * abs(e)) { // e is the second last p/2q
			absP = p;
			_DOUBLE_HI(absP) &= 0x7FFFFFFF;
			absQ = q * r;
			_DOUBLE_HI(absQ) &= 0x7FFFFFFF;

			if (absP < absQ) {
				// parabolic interpolation step
				d = p / (q + q);
				u = x + d;
				if (u < a || u > b) {
					// golden section interpolation step
					e = (x < m ? b : a) - x;
					d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
					////
					////count iters
					//golden_sec_iters++;
				}
				// f must not be evaluated too close to a or b
				else if ((u - a) < tol2 || (b - u) < tol2) {
					d = x < m
						? tol
						: -tol;
					////
					////count iters
					//parabolic_iters++;
				}
				//else
				//	//count iters
				//	parabolic_iters++;
			}
			else {
				// golden section interpolation step
				e = (x < m ? b : a) - x;
				d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
				////
				////count iters
				//golden_sec_iters++;
			}
		}
		else {
			// golden section interpolation step
			e = (x < m ? b : a) - x;
			d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
			////
			////count iters
			//golden_sec_iters++;
		}

		// f must not be evaluated too close to x
		absTemp = d;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp >= tol)
			u = x + d;
		else if (d > 0)
			u = x + tol;
		else
			u = x - tol;
		fu = func(u);

		// update Variables
		// save minimum
		if (fu <= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || x == w) {
				//update w and v if u lower than w
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || x == v || w == v) {
				//update v only if u lower than v
				v = u;
				fv = fu;
			}
		}

		--count;

		//update paramters
		m = a + (b - a) * 0.5;
		absTemp = x;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
		tol2 = tol + tol;

		absTemp = x - m;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	}

	//total_iters = golden_sec_iters + parabolic_iters;
	max_iter -= count;

	return new double[2] {x, fx}; //return x;
}

double* BrentMethodExtrema::FindMaxima(double (* const func)(double), double const left, double const right, unsigned int& max_iter) {
	//unsigned int iters = 0;
	//unsigned int total_iters = 0;
	//unsigned int golden_sec_iters = 0;
	//unsigned int parabolic_iters = 0;

	if (max_iter == 0) // special case
		return new double[2] {left, func(left)}; //return left;

	unsigned int count = max_iter;

	double a, b, m;
	if (left < right) {
		a = left;
		b = right;
	}
	else {
		a = right;
		b = left;
	}

	double x, w, v, u, fx, fw, fv, fu;

	v = w = x = a + (b - a) * _BRENTMETHOD_GOLDEN_RATIO_B;
	fv = fw = fx = func(x);

	double d, e;
	d = e = 0;

	double r, p, q;
	double absP, absQ;
	double absTemp;

	double tol, tol2;
	// distance between x and u, check tol
	// distance between x and boundaries, check tol2

	//update paramters
	m = a + (b - a) * 0.5;
	absTemp = x;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
	tol2 = tol + tol;

	absTemp = x - m;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;

	// f must not be evaluated too close to a or b
	while (count && absTemp > (tol2 - (b - a) * 0.5)) {
		//iters++;
		//check SPI behaving
		// p = (w - x) * (w - x) * (fx - fv) + (v - x) * (v - x) * (fw - fx);
		// q = (w - x) * (fx - fv) + (v - x) * (fw - fx);
		// u = x + p / (q + q);
		absTemp = e;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp > tol) {
			// fit parabolic
			r = (v - x) * (fw - fx);
			q = (w - x) * (fx - fv);
			p = (w - x) * q;
			q += r;
			p += (v - x) * r;
			//q += q;
			r = e; // r is the second last p/2q
			e = d; // e becomes last p/2q

			// if(abs(p/(2*q)) < 0.5 * abs(e)) { // e is the second last p/2q
			absP = p;
			_DOUBLE_HI(absP) &= 0x7FFFFFFF;
			absQ = q * r;
			_DOUBLE_HI(absQ) &= 0x7FFFFFFF;

			if (absP < absQ) {
				// parabolic interpolation step
				d = p / (q + q);
				u = x + d;
				if (u < a || u > b) {
					// golden section interpolation step
					e = (x < m ? b : a) - x;
					d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
					////
					////count iters
					//golden_sec_iters++;
				}
				// f must not be evaluated too close to a or b
				else if ((u - a) < tol2 || (b - u) < tol2) {
					d = x < m
						? tol
						: -tol;
					////
					////count iters
					//parabolic_iters++;
				}
				//else
				//	//count iters
				//	parabolic_iters++;
			}
			else {
				// golden section interpolation step
				e = (x < m ? b : a) - x;
				d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
				////
				////count iters
				//golden_sec_iters++;
			}
		}
		else {
			// golden section interpolation step
			e = (x < m ? b : a) - x;
			d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
			////
			////count iters
			//golden_sec_iters++;
		}

		// f must not be evaluated too close to x
		absTemp = d;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp >= tol)
			u = x + d;
		else if (d > 0)
			u = x + tol;
		else
			u = x - tol;
		fu = func(u);

		// update Variables
		// save maximum
		if (fu >= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu >= fw || x == w) {
				//update w and v if u greater than w
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu >= fv || x == v || w == v) {
				//update v only if u greater than v
				v = u;
				fv = fu;
			}
		}

		--count;

		//update paramters
		m = a + (b - a) * 0.5;
		absTemp = x;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
		tol2 = tol + tol;

		absTemp = x - m;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	}

	//total_iters = golden_sec_iters + parabolic_iters;
	max_iter -= count;

	return new double[2] {x, fx}; //return x;
}

double* BrentMethodExtrema::FindMinima(DoubleFunc* const func, double const left, double const right, unsigned int& max_iter) {
	//unsigned int iters = 0;
	//unsigned int total_iters = 0;
	//unsigned int golden_sec_iters = 0;
	//unsigned int parabolic_iters = 0;

	/**
	* x is the least value of function
	* w is the second lowest value of function
	* v is previous value of w. usually it is next the lowest value
	*
	* u is last evaluated point
	* a is left bound
	* b is right bound
	* m is interval midpoint
	*
	* tol is tolerance
	* tol2 is tolerance2
	*
	* -----------------------------------------------------------------/
	** /--parabolic interpolation variables--/
	* p is numerator
	* q is denominator
	* u = x + p/q
	*
	* d is distance between x and u, also equal p/q. The distance moved in the last step
	* e is the second last p/q. The distance moved in the step before last
	*
	** /--conditions--/
	* q != 0;
	* u belong to [a,b]
	* abs(p/q) < 0.5 * abs(e)
	* pass tolerance check
	*
	* condition "p/q < n" can change to "p < n*q", when q > 0.
	* changing to "p < n*q" can also avoid q == 0
	*/

	if (max_iter == 0) // special case
		return new double[2] {left, (*func)(left)}; //return left;

	unsigned int count = max_iter;

	double a, b, m;
	if (left < right) {
		a = left;
		b = right;
	}
	else {
		a = right;
		b = left;
	}

	double x, w, v, u, fx, fw, fv, fu;

	v = w = x = a + (b - a) * _BRENTMETHOD_GOLDEN_RATIO_B;
	fv = fw = fx = (*func)(x);

	double d, e;
	d = e = 0;

	double r, p, q;
	double absP, absQ;
	double absTemp;

	double tol, tol2;
	// distance between x and u, check tol
	// distance between x and boundaries, check tol2

	//update paramters
	m = a + (b - a) * 0.5;
	absTemp = x;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
	tol2 = tol + tol;

	absTemp = x - m;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;

	// f must not be evaluated too close to a or b
	while (count && absTemp > (tol2 - (b - a) * 0.5)) {
		//iters++;
		//check SPI behaving
		// p = (w - x) * (w - x) * (fx - fv) + (v - x) * (v - x) * (fw - fx);
		// q = (w - x) * (fx - fv) + (v - x) * (fw - fx);
		// u = x + p / (q + q);
		absTemp = e;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp > tol) {
			// fit parabolic
			r = (v - x) * (fw - fx);
			q = (w - x) * (fx - fv);
			p = (w - x) * q;
			q += r;
			p += (v - x) * r;
			//q += q;
			r = e; // r is the second last p/2q
			e = d; // e becomes last p/2q

			// if(abs(p/(2*q)) < 0.5 * abs(e)) { // e is the second last p/2q
			absP = p;
			_DOUBLE_HI(absP) &= 0x7FFFFFFF;
			absQ = q * r;
			_DOUBLE_HI(absQ) &= 0x7FFFFFFF;

			if (absP < absQ) {
				// parabolic interpolation step
				d = p / (q + q);
				u = x + d;
				if (u < a || u > b) {
					// golden section interpolation step
					e = (x < m ? b : a) - x;
					d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
					////
					////count iters
					//golden_sec_iters++;
				}
				// f must not be evaluated too close to a or b
				else if ((u - a) < tol2 || (b - u) < tol2) {
					d = x < m
						? tol
						: -tol;
					////
					////count iters
					//parabolic_iters++;
				}
				//else
				//	//count iters
				//	parabolic_iters++;
			}
			else {
				// golden section interpolation step
				e = (x < m ? b : a) - x;
				d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
				////
				////count iters
				//golden_sec_iters++;
			}
		}
		else {
			// golden section interpolation step
			e = (x < m ? b : a) - x;
			d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
			////
			////count iters
			//golden_sec_iters++;
		}

		// f must not be evaluated too close to x
		absTemp = d;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp >= tol)
			u = x + d;
		else if (d > 0)
			u = x + tol;
		else
			u = x - tol;
		fu = (*func)(u);

		// update Variables
		// save minimum
		if (fu <= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || x == w) {
				//update w and v if u lower than w
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || x == v || w == v) {
				//update v only if u lower than v
				v = u;
				fv = fu;
			}
		}

		--count;

		//update paramters
		m = a + (b - a) * 0.5;
		absTemp = x;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
		tol2 = tol + tol;

		absTemp = x - m;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	}

	//total_iters = golden_sec_iters + parabolic_iters;
	max_iter -= count;

	return new double[2] {x, fx}; //return x;
}

double* BrentMethodExtrema::FindMaxima(DoubleFunc* const func, double const left, double const right, unsigned int& max_iter) {
	//unsigned int iters = 0;
	//unsigned int total_iters = 0;
	//unsigned int golden_sec_iters = 0;
	//unsigned int parabolic_iters = 0;

	if (max_iter == 0) // special case
		return new double[2] {left, (*func)(left)}; //return left;

	unsigned int count = max_iter;

	double a, b, m;
	if (left < right) {
		a = left;
		b = right;
	}
	else {
		a = right;
		b = left;
	}

	double x, w, v, u, fx, fw, fv, fu;

	v = w = x = a + (b - a) * _BRENTMETHOD_GOLDEN_RATIO_B;
	fv = fw = fx = (*func)(x);

	double d, e;
	d = e = 0;

	double r, p, q;
	double absP, absQ;
	double absTemp;

	double tol, tol2;
	// distance between x and u, check tol
	// distance between x and boundaries, check tol2

	//update paramters
	m = a + (b - a) * 0.5;
	absTemp = x;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
	tol2 = tol + tol;

	absTemp = x - m;
	_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;

	// f must not be evaluated too close to a or b
	while (count && absTemp > (tol2 - (b - a) * 0.5)) {
		//iters++;
		//check SPI behaving
		// p = (w - x) * (w - x) * (fx - fv) + (v - x) * (v - x) * (fw - fx);
		// q = (w - x) * (fx - fv) + (v - x) * (fw - fx);
		// u = x + p / (q + q);
		absTemp = e;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp > tol) {
			// fit parabolic
			r = (v - x) * (fw - fx);
			q = (w - x) * (fx - fv);
			p = (w - x) * q;
			q += r;
			p += (v - x) * r;
			//q += q;
			r = e; // r is the second last p/2q
			e = d; // e becomes last p/2q

			// if(abs(p/(2*q)) < 0.5 * abs(e)) { // e is the second last p/2q
			absP = p;
			_DOUBLE_HI(absP) &= 0x7FFFFFFF;
			absQ = q * r;
			_DOUBLE_HI(absQ) &= 0x7FFFFFFF;

			if (absP < absQ) {
				// parabolic interpolation step
				d = p / (q + q);
				u = x + d;
				if (u < a || u > b) {
					// golden section interpolation step
					e = (x < m ? b : a) - x;
					d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
					////
					////count iters
					//golden_sec_iters++;
				}
				// f must not be evaluated too close to a or b
				else if ((u - a) < tol2 || (b - u) < tol2) {
					d = x < m
						? tol
						: -tol;
					////
					////count iters
					//parabolic_iters++;
				}
				//else
				//	//count iters
				//	parabolic_iters++;
			}
			else {
				// golden section interpolation step
				e = (x < m ? b : a) - x;
				d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
				////
				////count iters
				//golden_sec_iters++;
			}
		}
		else {
			// golden section interpolation step
			e = (x < m ? b : a) - x;
			d = e * _BRENTMETHOD_GOLDEN_RATIO_B;
			////
			////count iters
			//golden_sec_iters++;
		}

		// f must not be evaluated too close to x
		absTemp = d;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		if (absTemp >= tol)
			u = x + d;
		else if (d > 0)
			u = x + tol;
		else
			u = x - tol;
		fu = (*func)(u);

		// update Variables
		// save maximum
		if (fu >= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu >= fw || x == w) {
				//update w and v if u greater than w
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu >= fv || x == v || w == v) {
				//update v only if u greater than v
				v = u;
				fv = fu;
			}
		}

		--count;

		//update paramters
		m = a + (b - a) * 0.5;
		absTemp = x;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
		tol = _BRENTMETHOD_SQRT_DOUBLE_EPSILON * absTemp + _BRENTMETHOD_SQRT_DOUBLE_EPSILON_QUARTER; // t = _SQRT_EPSILON_QUARTER, t must be more than 0 to avoid 0 tolerance
		tol2 = tol + tol;

		absTemp = x - m;
		_DOUBLE_HI(absTemp) &= 0x7FFFFFFF;
	}

	//total_iters = golden_sec_iters + parabolic_iters;
	max_iter -= count;

	return new double[2] {x, fx}; //return x;
}

#include "Remez.h"

using namespace MyMath;

/****
* references:
* boost.org, .\libs\math\include_private\boost\math\tools\remez.hpp
*/

Remez::Remez(double (* const func)(double)
    , const unsigned int oN, const unsigned int oD
    , const double left_interval, const double right_interval
    , const int skew)
    :Remez(func
        , true, false
        , oN, oD
        , left_interval, right_interval
        , skew) {
}

Remez::Remez(double (* const func)(double)
    , const bool pinned
    , const unsigned int oN, const unsigned int oD
    , const double left_interval, const double right_interval
    , const int skew)
	:Remez(func
		, true, pinned
		, oN, oD
		, left_interval, right_interval
		, skew) {
}

Remez::Remez(double (* const func)(double)
    , const bool relative_error, const bool pinned
    , const unsigned int oN, const unsigned int oD
    , const double left_interval, const double right_interval
    , const int skew) {
    //initialize
    Remez::func = func;
    Remez::weight_func = NULL;

    //Remez::relative_error = relative_error;
    Remez::error_type = relative_error
        ? RemezErrorType::RelativeError
        : RemezErrorType::AbsoluteError;
    Remez::pinned = pinned;

    Remez::oN = oN;
    Remez::oD = oD;

    if (left_interval < right_interval) {
        Remez::left_interval = left_interval, Remez::right_interval = right_interval;
    }
    else {
        Remez::left_interval = right_interval, Remez::right_interval = left_interval;
    }

    Remez::skew = skew;
    Remez::brake = 0;

    maxRank = oN > oD ? oN : oD;

    n = pinned ? (oN + oD) : (oN + oD + 1);
    unknownCount = n + 1; // add 1 for Error at remez step
    rootCount = n + 2;

    Remez::iter_count = 0;

    Remez::sanity = true;

    Remez::solution = new double[unknownCount];
    Remez::matrix = new double* [unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        Remez::matrix[i] = new double[unknownCount];

    Remez::estimFuncC = new EstimateFuncContainer{ solution , pinned , pinned ? oN : (oN + 1) , oD };
    Remez::funcC = new FuncContainer{ func };
    Remez::weight_funcC = NULL;

    Remez::error_func_roots = Remez::GenerateChebyshevKnotsAndBoundaries();// descending order
    Remez::error_func_extrema = new double* [unknownCount]; // add 1 for Error at remez step

    Remez::max_error = 0.0;
    Remez::max_error_change = 0.0;
    Remez::max_error_last_change = 0.0;

    Remez::solution_max_abs_error = 0.0;
    Remez::solution_max_rel_error = 0.0;

    Remez::Initialize();
}

Remez::Remez(double (*func)(double)
    , double (*weight_func)(double), bool pinned
    , unsigned int oN, unsigned int oD
    , double left_interval, double right_interval
    , int skew) {
    //initialize
    Remez::func = func;
    Remez::weight_func = weight_func;

    //Remez::relative_error = false;
    Remez::error_type = weight_func
        ? RemezErrorType::WeightError
        : RemezErrorType::AbsoluteError;
    Remez::pinned = pinned;

    Remez::oN = oN;
    Remez::oD = oD;

    if (left_interval < right_interval) {
        Remez::left_interval = left_interval, Remez::right_interval = right_interval;
    }
    else {
        Remez::left_interval = right_interval, Remez::right_interval = left_interval;
    }

    Remez::skew = skew;
    Remez::brake = 0;

    maxRank = oN > oD ? oN : oD;

    n = pinned ? (oN + oD) : (oN + oD + 1);
    unknownCount = n + 1; // add 1 for Error at remez step
    rootCount = n + 2;

    Remez::iter_count = 0;

    Remez::sanity = true;

    Remez::solution = new double[unknownCount];
    Remez::matrix = new double* [unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        Remez::matrix[i] = new double[unknownCount];

    Remez::estimFuncC = new EstimateFuncContainer{ solution , pinned , pinned ? oN : (oN + 1) , oD };
    Remez::funcC = new FuncContainer{ func };
    if (Remez::error_type == RemezErrorType::WeightError)
        Remez::weight_funcC = new FuncContainer{ weight_func };
    else
        Remez::weight_funcC = NULL;

    Remez::error_func_roots = Remez::GenerateChebyshevKnotsAndBoundaries();// descending order
    Remez::error_func_extrema = new double* [unknownCount]; // add 1 for Error at remez step

    Remez::max_error = 0.0;
    Remez::max_error_change = 0.0;
    Remez::max_error_last_change = 0.0;

    Remez::solution_max_abs_error = 0.0;
    Remez::solution_max_rel_error = 0.0;

    Remez::Initialize();
}

Remez::~Remez() {
	delete[] solution;
	for (unsigned int i = 0; i < unknownCount; i++)
		delete[] matrix[i];
    delete[] matrix;

    delete estimFuncC;
    delete funcC;
    delete weight_funcC;

	delete[] error_func_roots;
    for (unsigned int i = 0; i < unknownCount; i++)
        delete[] error_func_extrema[i];
    delete[] error_func_extrema;
}

unsigned int Remez::GetBrake() {
    return Remez::brake;
}

void Remez::SetBrake(unsigned int brake) {
	Remez::brake = brake;
}

bool Remez::CheckValidNum(double& x) {
    return (_DOUBLE_HI(x) & 0x7FF00000) != 0x7FF00000; // check if exp bits are all 1, then it may be +inf, -inf or +-NaN.
}

double* Remez::GenerateChebyshevKnotsAndBoundaries() { // return double array
	/**
	* first kind
	* 1. original Chebyshev knots in interval [-1, 1] :
	* x = cos( (2*k + 1) * pi/(2*n) ); where k = 0,...,n
	* 
	* 2. change interval to [0, 1]
	* x = 0.5 * (cos( (2*k + 1) * pi/(2*n) ) + 1); where k = 0,...,n
	* 
	* 3. skew, positive for left, negative for right
	* x = (0.5 * (cos( (2*k + 1) * pi/(2*n) ) + 1))**(skew/200 + 1); where k = 0,...,n; skew in [-100, 100]
	* 
	* 4. rescale over interval [a, b], a = left bound, b = right bound
	* x = a + (b - a) * (0.5 * (cos( (2*k + 1) * pi/(2*n) ) + 1))**(skew/200 + 1); where k = 0,...,n; skew in [-100, 100], a < b
	* f(x) = cos(n*arccos( 2*( (x-a)/(b-a) )**( 200/(200 + skew) ) - 1 ))
    * 
    * demo URL:
    * https://www.desmos.com/calculator/h4v2ah7uvq
    * 
    * n = count of roots
    * i = index of x
    * a = left bound
    * b = right bound
    * s = skew
	*/

	if (rootCount < 1)
		return NULL;

	double* x = new double[rootCount];
	//unsigned int n;
	double c, lenght, k;
	if (skew == 0) {
		switch (rootCount) {
		default:
			//n = count - 2;
			c = PI / (n + n);
			lenght = right_interval - left_interval;
			k = 0.0;
			for (unsigned int i = 1; i <= n; i++, k++)
				x[i] = (cos((k + k + 1) * c) + 1) * 0.5 * lenght + left_interval;
		case 2:
			x[rootCount - 1] = left_interval;
		case 1:
			x[0] = right_interval;
		}
	}
	else {
		double skewRatioPower; // skew/200 + 1
		switch (rootCount) {
		default:
			//n = count - 2;
			c = PI / (n + n);
			lenght = right_interval - left_interval;
			k = 0.0;
			skewRatioPower = skew * 0.005 + 1; // skew/200 + 1
			for (unsigned int i = 1; i <= n; i++, k++)
				x[i] = pow((cos((k + k + 1) * c) + 1) * 0.5, skewRatioPower) * lenght + left_interval;
		case 2:
			x[rootCount - 1] = left_interval;
		case 1:
			x[0] = right_interval;
		}
	}

	return x;
}

void Remez::Initialize() {

	//do first approximate
	// 
	// 1. P(x) / [1 + Q(x)] = f(x)
	// 2. P(x) = f(x)*[1 + Q(x)]
	// 3. P(x) - f(x)*Q(x) = f(x)
	// 
	// Cm for numerator's coefficients, Cn for denominator's coefficients
	// f(x) = Cm0 + Cm1 * x + Cm2 * x**2 + ... + CmoN * x**oN
	//  - f(x) * (Cn1 * x + Cn2 * x**2 + ... + CnoD * x**oD)
	// all Cm and Cn are unknowns
	// if pinned then Cm0 is 0, else Cm0 is unknown
	// b = Ax
	//build b vector

	for (unsigned int i = 0, i2 = 1; i < n; i++, i2++) {
		solution[i] = func(error_func_roots[i2]); // function
	}
	solution[n] = 0; // insert 0 error for guess E step

	//build n x n matrix
	if (pinned) {
		double x, xx;
		for (unsigned int i = 0, i2 = 1; i < n; i++, i2++) {
			xx = x = error_func_roots[i2];
			for (unsigned int j = 0; j < maxRank; j++) {
				if (j < oN)
					matrix[i][j] = xx; // numerator
				if (j + oN < n)
					matrix[i][j + oN] = -xx * solution[i]; //-xx * func(x); // denominator

				xx *= x;
			}
		}
	}
	else {
		double x, xx;
		for (unsigned int i = 0, i2 = 1; i < n; i++, i2++) {
			xx = x = error_func_roots[i2];
			matrix[i][0] = 1; // Cm0
			for (unsigned int j = 1; j <= maxRank; j++) {
				if (j <= oN)
					matrix[i][j] = xx; // numerator
				if (j + oN < n)
					matrix[i][j + oN] = -xx * solution[i]; //-xx * func(x); // denominator

				xx *= x;
			}
		}
	}

    LUP::SolvingLinearEquations(n, matrix, solution); // n x n

    //check whether solution is valid
    for (unsigned int i = 0; i < unknownCount; i++) {
        //check whether value is valid
        if (!CheckValidNum(solution[i])) {
            sanity = false;
            break;
        }
    }

	// update extrema
	unsigned int iter;
    switch (error_type) {
    case RemezErrorType::WeightError:
        InitializeWeightedErrorFuncExtrema(iter);
        break;
    case RemezErrorType::RelativeError:
        InitializeRelativeErrorFuncExtrema(iter);
        break;
    default: //case RemezErrorType::AbsoluteError:
        InitializeAbsoluteErrorFuncExtrema(iter);
        break;
    }

    max_error_change = max_error;
}

bool Remez::Iterate() { // return success
    if (!sanity)
        return false;

    //iterate

    double E = solution[n];
    double* bVector = new double[unknownCount];
    double* yVector = new double[unknownCount];

    switch (error_type) {
    case RemezErrorType::WeightError:
        IteratingSolveLinearEquationsWithWeightError(E, bVector, yVector);
        break;
    case RemezErrorType::RelativeError:
        IteratingSolveLinearEquationsWithRelativeError(E, bVector, yVector);
        break;
    case RemezErrorType::AbsoluteError:
        IteratingSolveLinearEquationsWithAbsoluteError(E, bVector, yVector);
        break;
    }

    iter_count++;

    //check whether solution is valid
    if (!CheckSolutionValid(bVector)) {
        // fail
        delete[] bVector;
        delete[] yVector;

        return sanity = false;
    }

    //check solution error
    if (!CheckSolutionError(bVector, yVector)) {
        // too much error
        delete[] bVector;
        delete[] yVector;

        return sanity = false;
    }

    //save solution
    for (unsigned int i = 0; i < unknownCount; i++)
        solution[i] = bVector[i];

    //update control points
    //check control points alternate in sign first, if not, try to correct it
    if (!CheckErrorFuncExtremaAlternateSign()) {
        //fatal failed, cannot correct it
        delete[] bVector;
        delete[] yVector;

        return sanity = false;
    }

    ////check control points alternate in sign agian
    //double e1, e2;
    //AbsoluteErrorFunc errorFunc = { estimFuncC, funcC };
    //for (unsigned int i = 1; i <= n; i++) {
    //    e1 = errorFunc(controlPoints[i - 1][0]);
    //    e2 = errorFunc(controlPoints[i][0]);
    //    // check sign
    //    if ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)) {
    //        for (unsigned int i = 0; i < unknownCount; i++)
    //            delete[] AMatrix[i];
    //        delete[] AMatrix;

    //        delete[] bVector;
    //        delete[] yVector;

    //        return sanity = false;
    //    }
    //}

    unsigned int iter;
    // update roots
    switch (error_type) {
    case RemezErrorType::WeightError:
        UpdateWeightedErrorFuncRoots(iter);
        break;
    case RemezErrorType::RelativeError:
        UpdateRelativeErrorFuncRoots(iter);
        break;
    case RemezErrorType::AbsoluteError:
        UpdateAbsoluteErrorFuncRoots(iter);
        break;
    }

    max_error_last_change = max_error_change;
    max_error_change = max_error;

    // update extrema
    if (brake) {
        // controlPoint = newExt_x + brake% * (oldExt_x - newExt_x)
        const double brakePercent = (double)brake / 100.0;
        switch (error_type) {
        case RemezErrorType::WeightError:
            UpdateWeightedErrorFuncExtrema(brakePercent, iter);
            break;
        case RemezErrorType::RelativeError:
            UpdateRelativeErrorFuncExtrema(brakePercent, iter);
            break;
        default: //case RemezErrorType::AbsoluteError:
            UpdateAbsoluteErrorFuncExtrema(brakePercent, iter);
            break;
        }
    }
    else { // brake == 0
        switch (error_type) {
        case RemezErrorType::WeightError:
            UpdateWeightedErrorFuncExtrema(iter);
            break;
        case RemezErrorType::RelativeError:
            UpdateRelativeErrorFuncExtrema(iter);
            break;
        default: //case RemezErrorType::AbsoluteError:
            UpdateAbsoluteErrorFuncExtrema(iter);
            break;
        }
    }

    max_error_change = max_error - max_error_change;

    delete[] bVector;
    delete[] yVector;

    return sanity = true;
}

void Remez::IteratingSolveLinearEquationsWithWeightError(double E, double*& solutionVector, double*& funcVector) {
    double**& controlPoints = error_func_extrema;
    double*& bVector = solutionVector;
    double*& yVector = funcVector;

    // f(x) - P(x) / [1 + Q(x)] = E
    // 1. P(x) / [1 + Q(x)] + E = f(x)
    // 2. P(x) + E*[1 + Q(x)] = f(x)*[1 + Q(x)]
    // 3. P(x) - f(x)*Q(x) + E*Q(x) + E = f(x)
    // 4. P(x) - [f(x) - E]*Q(x) + E = f(x)
    // 5. P(x) + [E - f(x)]*Q(x) + E = f(x)
    // **first iterate
    // current truncated polynomial in first iterate is "P(x) + [E - f(x)]*Q(x)"
    // **second iterate
    // because E is also an unknown, need to approximate E first
    // so the equation in second iterate is "P(x) + [E' - f(x)]*Q(x) + E = f(x)", E' is the last E in last iterate
    // **need to guess real E first
    // 
    // * for using weight function
    // {f(x) - P(x) / [1 + Q(x)]} / W(x) = E
    // f(x) - P(x) / [1 + Q(x)] = E * W(x)
    // ...
    // P(x) + [(E' * W(x)) - f(x)]*Q(x) + (E * W(x)) = f(x)
    // 
    // Cm for numerator's coefficients, Cn for denominator's coefficients
    // f(x) = Cm0 + Cm1 * x + Cm2 * x**2 + ... + CmoN * x**oN
    //  + [E' - f(x)] * (Cn1 * x + Cn2 * x**2 + ... + CnoD * x**oD)
    //  + E
    // all Cm, Cn and E are unknowns
    // if pinned then Cm0 is 0, else Cm0 is unknown
    // 

    double sign;
    double guessE;
    double err_err; // measure error's error
    double err_err_tol = _REMEZ_DEFAULT_RATIONAL_ERR_TOLERANCE; // error's error tolerance
    unsigned int err_max_convergence = _REMEZ_DEFAULT_RATIONAL_ERR_MAX_CONVGC;

    double** AMatrix = new double* [unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        AMatrix[i] = new double[unknownCount];

    double x, xx, y;
    double absGuessE, absE = E;
    _DOUBLE_HI(absE) &= 0x7FFFFFFF;

    double w;
    //pre-calculate weights
    double* wVector = new double[unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        wVector[i] = weight_func(controlPoints[i][0]); // function

    if (pinned) {
        // if "pinned" is true, all control points cannot be 0. The origin(0,0) control point will occur the 0 dividing during LUP step.
        for (unsigned int i = 0; i < unknownCount; i++) {
            //avoid the origin(0,0) control point
            if (controlPoints[i][0] == 0.0) {
                controlPoints[i][0] = controlPoints[i == 0 ? (i + 1) : (i - 1)][0] * _BRENTMETHOD_GOLDEN_RATIO_B; // /3
                controlPoints[i][1] = yVector[i] = func(controlPoints[i][0]);
                _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
            }
            else
                yVector[i] = func(controlPoints[i][0]); // function
        }
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                y = bVector[i] = yVector[i];
                w = wVector[i];
                for (unsigned int j = 0; j < maxRank; j++) {
                    if (j < oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E * w - y); //-xx * ( func(x) - E * weight_func(x) ); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign * w;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }
    else {
        for (unsigned int i = 0; i < unknownCount; i++)
            yVector[i] = func(controlPoints[i][0]); // function
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                y = bVector[i] = yVector[i];
                w = wVector[i];
                matrix[i][0] = AMatrix[i][0] = 1; // Cm0
                for (unsigned int j = 1; j <= maxRank; j++) {
                    if (j <= oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E * w - y); //-xx * ( func(x) - E * weight_func(x) ); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign * w;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }

    delete[] wVector;
    for (unsigned int i = 0; i < unknownCount; i++)
        delete[] AMatrix[i];
    delete[] AMatrix;
}

void Remez::IteratingSolveLinearEquationsWithRelativeError(double E, double*& solutionVector, double*& funcVector) {
    double**& controlPoints = error_func_extrema;
    double*& bVector = solutionVector;
    double*& yVector = funcVector;

    // f(x) - P(x) / [1 + Q(x)] = E
    // 1. P(x) / [1 + Q(x)] + E = f(x)
    // 2. P(x) + E*[1 + Q(x)] = f(x)*[1 + Q(x)]
    // 3. P(x) - f(x)*Q(x) + E*Q(x) + E = f(x)
    // 4. P(x) - [f(x) - E]*Q(x) + E = f(x)
    // 5. P(x) + [E - f(x)]*Q(x) + E = f(x)
    // **first iterate
    // current truncated polynomial in first iterate is "P(x) + [E - f(x)]*Q(x)"
    // **second iterate
    // because E is also an unknown, need to approximate E first
    // so the equation in second iterate is "P(x) + [E' - f(x)]*Q(x) + E = f(x)", E' is the last E in last iterate
    // **need to guess real E first
    // 
    // * for using relative error
    // {f(x) - P(x) / [1 + Q(x)]} / |f(x)| = E
    // f(x) - P(x) / [1 + Q(x)] = E * |f(x)|
    // ...
    // P(x) + [(E' * |f(x)|) - f(x)]*Q(x) + (E * |f(x)|) = f(x)
    // 
    // Cm for numerator's coefficients, Cn for denominator's coefficients
    // f(x) = Cm0 + Cm1 * x + Cm2 * x**2 + ... + CmoN * x**oN
    //  + [E' - f(x)] * (Cn1 * x + Cn2 * x**2 + ... + CnoD * x**oD)
    //  + E
    // all Cm, Cn and E are unknowns
    // if pinned then Cm0 is 0, else Cm0 is unknown
    // 

    double sign;
    double guessE;
    double err_err; // measure error's error
    double err_err_tol = _REMEZ_DEFAULT_RATIONAL_ERR_TOLERANCE; // error's error tolerance
    unsigned int err_max_convergence = _REMEZ_DEFAULT_RATIONAL_ERR_MAX_CONVGC;

    double** AMatrix = new double* [unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        AMatrix[i] = new double[unknownCount];

    double x, xx, y;
    double absGuessE, absE = E;
    _DOUBLE_HI(absE) &= 0x7FFFFFFF;

    double absY;

    if (pinned) {
        // if "pinned" is true, all control points cannot be 0. The origin(0,0) control point will occur the 0 dividing during LUP step.
        for (unsigned int i = 0; i < unknownCount; i++) {
            //avoid the origin(0,0) control point
            if (controlPoints[i][0] == 0.0) {
                controlPoints[i][0] = controlPoints[i == 0 ? (i + 1) : (i - 1)][0] * _BRENTMETHOD_GOLDEN_RATIO_B; // /3
                controlPoints[i][1] = yVector[i] = func(controlPoints[i][0]);
                _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
            }
            else
                yVector[i] = func(controlPoints[i][0]); // function
        }
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                absY = y = bVector[i] = yVector[i];
                _DOUBLE_HI(absY) &= 0x7FFFFFFF;
                for (unsigned int j = 0; j < maxRank; j++) {
                    if (j < oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E * absY - y); //-xx * ( func(x) - E * abs(func(x)) ); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign * absY;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }
    else {
        for (unsigned int i = 0; i < unknownCount; i++)
            yVector[i] = func(controlPoints[i][0]); // function
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                absY = y = bVector[i] = yVector[i];
                _DOUBLE_HI(absY) &= 0x7FFFFFFF;
                matrix[i][0] = AMatrix[i][0] = 1; // Cm0
                for (unsigned int j = 1; j <= maxRank; j++) {
                    if (j <= oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E * absY - y); //-xx * ( func(x) - E * abs(func(x)) ); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign * absY;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }

    for (unsigned int i = 0; i < unknownCount; i++)
        delete[] AMatrix[i];
    delete[] AMatrix;
}

void Remez::IteratingSolveLinearEquationsWithAbsoluteError(double E, double*& solutionVector, double*& funcVector) {
    double**& controlPoints = error_func_extrema;
    double*& bVector = solutionVector;
    double*& yVector = funcVector;
    //
    double sign;
    double guessE;
    double err_err; // measure error's error
    double err_err_tol = _REMEZ_DEFAULT_RATIONAL_ERR_TOLERANCE; // error's error tolerance
    unsigned int err_max_convergence = _REMEZ_DEFAULT_RATIONAL_ERR_MAX_CONVGC;

    double** AMatrix = new double* [unknownCount];
    for (unsigned int i = 0; i < unknownCount; i++)
        AMatrix[i] = new double[unknownCount];

    double x, xx, y;
    double absGuessE, absE = E;
    _DOUBLE_HI(absE) &= 0x7FFFFFFF;

    if (pinned) {
        // if "pinned" is true, all control points cannot be 0. The origin(0,0) control point will occur the 0 dividing during LUP step.
        for (unsigned int i = 0; i < unknownCount; i++) {
            //avoid the origin(0,0) control point
            if (controlPoints[i][0] == 0.0) {
                controlPoints[i][0] = controlPoints[i == 0 ? (i + 1) : (i - 1)][0] * _BRENTMETHOD_GOLDEN_RATIO_B; // /3
                controlPoints[i][1] = yVector[i] = func(controlPoints[i][0]);
                _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
            }
            else
                yVector[i] = func(controlPoints[i][0]); // function
        }
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                y = bVector[i] = yVector[i];
                for (unsigned int j = 0; j < maxRank; j++) {
                    if (j < oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E - y); //-xx * (func(x) - E); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }
    else {
        for (unsigned int i = 0; i < unknownCount; i++)
            yVector[i] = func(controlPoints[i][0]); // function
        //guess E
        do {
            sign = 1;
            for (unsigned int i = 0; i < unknownCount; i++) {
                xx = x = controlPoints[i][0];
                y = bVector[i] = yVector[i];
                matrix[i][0] = AMatrix[i][0] = 1; // Cm0
                for (unsigned int j = 1; j <= maxRank; j++) {
                    if (j <= oN)
                        matrix[i][j] = AMatrix[i][j] = xx; // numerator
                    if (j + oN < n)
                        matrix[i][j + oN] = AMatrix[i][j + oN] = xx * (E - y); //-xx * (func(x) - E); // denominator

                    xx *= x;
                }
                matrix[i][n] = AMatrix[i][n] = sign;
                _DOUBLE_HI(sign) ^= 0x80000000; // change sign bit
                _DOUBLE_HI(E) ^= 0x80000000; // change sign bit
            }

            LUP::SolvingLinearEquations(unknownCount, AMatrix, bVector);

            absGuessE = guessE = bVector[n];
            _DOUBLE_HI(absGuessE) &= 0x7FFFFFFF;

            if (E == 0) {
                err_err = 1.0;
            }
            else {
                err_err = (absGuessE - absE) / absE;
                _DOUBLE_HI(err_err) &= 0x7FFFFFFF;
            }

            absE = E = guessE;
            _DOUBLE_HI(absE) &= 0x7FFFFFFF;

        } while (oD // check if it is approximating by a rational function or a polynomial
            && --err_max_convergence
            && (err_err > err_err_tol));
    }

    for (unsigned int i = 0; i < unknownCount; i++)
        delete[] AMatrix[i];
    delete[] AMatrix;
}

//check sanity
bool Remez::CheckSolutionValid(double*& solutionVector) {
    for (unsigned int i = 0; i < unknownCount; i++) {
        //check whether value is valid
        if (!CheckValidNum(solutionVector[i]))
            return false;
    }
    return true;
}

bool Remez::CheckSolutionError(double*& solutionVector, double*& funcVector) {
    //Ax = y
    solution_max_abs_error = 0.0;
    solution_max_rel_error = 0.0;
    double solution_abs_error;
    double solution_rel_error;
    for (unsigned int i = 0; i < unknownCount; i++) {
        solution_abs_error = 0.0;
        for (unsigned int j = 0; j < unknownCount; j++) {
            solution_abs_error += matrix[i][j] * solutionVector[j];
        }
        solution_abs_error = solution_abs_error - funcVector[i];
        _DOUBLE_HI(solution_abs_error) &= 0x7FFFFFFF;

        solution_rel_error = funcVector[i];
        _DOUBLE_HI(solution_rel_error) &= 0x7FFFFFFF;
        solution_rel_error = solution_abs_error / solution_rel_error;

        if (solution_abs_error > solution_max_abs_error)
            solution_max_abs_error = solution_abs_error; // near 0 compare
        if (solution_rel_error > solution_max_rel_error)
            solution_max_rel_error = solution_rel_error; // large numbers compare
    }

    return solution_max_rel_error <= _REMEZ_SQRT_DOUBLE_EPSILON;
}

bool Remez::CheckErrorFuncExtremaAlternateSign() {
    double**& controlPoints = error_func_extrema;
    //check control points alternate in sign first
    double e1, e2, newPoint_x, perturbation;
    AbsoluteErrorFunc errorFunc = { estimFuncC, funcC };
    //
    e1 = errorFunc(controlPoints[0][0]);
    for (unsigned int i = 1; i < n; i++) {
        //e1 = errorFunc(controlPoints[i - 1][0]);
        e2 = errorFunc(controlPoints[i][0]);
        // check sign
        if ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)) {
            //try to find other different sign point to correct this failure
            /*
            * demo URL :
            * https://www.desmos.com/calculator/kb3sr4nk68
            *
            * c = current control point
            * L = left control point
            * R = right control point
            * p = perturbation
            */
            perturbation = _REMEZ_DEFAULT_CTRLPT_PERTURBATION;
            // current to left
            do {
                // newPoint_x = oldCurrentPoint_x + perturbation * (oldLeftPoint_x - oldCurrentPoint_x)
                newPoint_x = controlPoints[i][0] + (controlPoints[i - 1][0] - controlPoints[i][0]) * perturbation;
                e2 = errorFunc(newPoint_x);
            } while ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)
                && (perturbation += _REMEZ_DEFAULT_CTRLPT_PERTURBATION) < _REMEZ_DEFAULT_CTRLPT_MAX_PERTURBATION);

            if ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)) {
                perturbation = _REMEZ_DEFAULT_CTRLPT_PERTURBATION;
                // current to right
                do {
                    // newPoint_x = oldCurrentPoint_x + perturbation * (oldRightPoint_x - oldCurrentPoint_x)
                    newPoint_x = controlPoints[i][0] + (controlPoints[i + 1][0] - controlPoints[i][0]) * perturbation;
                    e2 = errorFunc(newPoint_x);
                } while ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)
                    && (perturbation += _REMEZ_DEFAULT_CTRLPT_PERTURBATION) < _REMEZ_DEFAULT_CTRLPT_MAX_PERTURBATION);
            }

            if ((_DOUBLE_HI(e1) & 0x80000000) == (_DOUBLE_HI(e2) & 0x80000000)) {
                //fatal failed, cannot correct it
                return false;
            }
            else {
                controlPoints[i][0] = newPoint_x;
                controlPoints[i][1] = e2;
                _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
            }
        }
        e1 = e2;
    }

    return true;
}

//update extrema or roots
void Remez::UpdateWeightedErrorFuncRoots(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    WeightedErrorFunc errFunc = { estimFuncC, funcC, weight_funcC };
    for (unsigned int i = 1; i <= n; i++) {
        iter = _REMEZ_DEFAULT_ROOTS_MAX_ITER;
        roots[i] = TOMS748Algorithm::FindRoot(&errFunc, controlPoints[i - 1][0], controlPoints[i][0], iter);
    }
}

void Remez::UpdateRelativeErrorFuncRoots(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    RelativeErrorFunc errFunc = { estimFuncC, funcC };
    for (unsigned int i = 1; i <= n; i++) {
        iter = _REMEZ_DEFAULT_ROOTS_MAX_ITER;
        roots[i] = TOMS748Algorithm::FindRoot(&errFunc, controlPoints[i - 1][0], controlPoints[i][0], iter);
    }
}

void Remez::UpdateAbsoluteErrorFuncRoots(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    AbsoluteErrorFunc errFunc = { estimFuncC, funcC };
    for (unsigned int i = 1; i <= n; i++) {
        iter = _REMEZ_DEFAULT_ROOTS_MAX_ITER;
        roots[i] = TOMS748Algorithm::FindRoot(&errFunc, controlPoints[i - 1][0], controlPoints[i][0], iter);
    }
}

void Remez::UpdateWeightedErrorFuncExtrema(const double brakePercent, unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_WeightedErrorFunc abs_errFunc = { estimFuncC, funcC, weight_funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;

        double* newCtrlPt = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        controlPoints[i][0] = newCtrlPt[0] + (controlPoints[i][0] - newCtrlPt[0]) * brakePercent;
        controlPoints[i][1] = func(controlPoints[i][0]);
        _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
        delete[] newCtrlPt;

        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::UpdateRelativeErrorFuncExtrema(const double brakePercent, unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_RelativeErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;

        double* newCtrlPt = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        controlPoints[i][0] = newCtrlPt[0] + (controlPoints[i][0] - newCtrlPt[0]) * brakePercent;
        controlPoints[i][1] = func(controlPoints[i][0]);
        _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
        delete[] newCtrlPt;

        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::UpdateAbsoluteErrorFuncExtrema(const double brakePercent, unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_AbsoluteErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;

        double* newCtrlPt = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        controlPoints[i][0] = newCtrlPt[0] + (controlPoints[i][0] - newCtrlPt[0]) * brakePercent;
        controlPoints[i][1] = func(controlPoints[i][0]);
        _DOUBLE_HI(controlPoints[i][1]) &= 0x7FFFFFFF;
        delete[] newCtrlPt;

        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::UpdateWeightedErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_WeightedErrorFunc abs_errFunc = { estimFuncC, funcC, weight_funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::UpdateRelativeErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_RelativeErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::UpdateAbsoluteErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_AbsoluteErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::InitializeWeightedErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_WeightedErrorFunc abs_errFunc = { estimFuncC, funcC, weight_funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::InitializeRelativeErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_RelativeErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

void Remez::InitializeAbsoluteErrorFuncExtrema(unsigned int& iter) {
    double*& roots = error_func_roots;
    double**& controlPoints = error_func_extrema;
    max_error = 0.0;
    Abs_AbsoluteErrorFunc abs_errFunc = { estimFuncC, funcC };
    for (unsigned int i = 0; i < unknownCount; i++) {
        //delete[] controlPoints[i];
        iter = _REMEZ_DEFAULT_EXTREMA_MAX_ITER;
        controlPoints[i] = BrentMethodExtrema::FindMaxima(&abs_errFunc, roots[i + 1], roots[i], iter);
        if (controlPoints[i][1] > max_error)
            max_error = controlPoints[i][1];
    }
}

RemezErrorType Remez::GetErrorType() {
    return Remez::error_type;
}

double Remez::GetMaxError() {
    return Remez::max_error;
}

double Remez::GetCurrentChangeOfMaxError() {
    return Remez::max_error_change;
}

double Remez::GetLastChangeOfMaxError() {
    return Remez::max_error_last_change;
}

double Remez::GetSolutionMaxAbsoluteError() {
    return Remez::solution_max_abs_error;
}

double Remez::GetSolutionMaxRelativeError() {
    return Remez::solution_max_rel_error;
}

unsigned int Remez::GetIterationCount() {
    return Remez::iter_count;
}

bool Remez::Sanity() {
    return Remez::sanity;
}

unsigned int Remez::GetErrorFuncRoots(double*& roots) {
    if (n > 0)
        roots = &(error_func_roots[1]);
    else
        roots = NULL;
    return n;
}

unsigned int Remez::GetControlPoints(double**& extrema) {
    if (n > 0)
        extrema = &(error_func_extrema[0]);
    else
        extrema = NULL;
    return n + 1;
}

unsigned int Remez::GetNumerator(double*& coefficients) {
    if (pinned) {
        if (oN > 0)
            coefficients = &(solution[0]);
        else
            coefficients = NULL;

        return oN;
    }
    else {
        //if ((oN + 1) > 0)
        coefficients = &(solution[0]);
        //else
        //    coefficients = NULL;

        return oN + 1;
    }
}

unsigned int Remez::GetDenominator(double*& coefficients) {
    if (oD > 0)
        coefficients = &(solution[pinned ? oN : (oN + 1)]);
    else
        coefficients = NULL;

    return oD;
}

double Remez::Estimate(double x) {
    return (*estimFuncC)(x);
}

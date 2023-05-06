#pragma once

#include <math.h>
#include "BrentMethodExtrema.h"
#include "TOMS748Algorithm.h"
#include "LUP.h"
#include "ErrorFunc.h"

#ifndef _DOUBLE_HI
#define _DOUBLE_HI(x) *(1 + (int*)&x)
#endif

#ifndef PI
#define PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

#ifndef _REMEZ_DEFAULT
#define _REMEZ_DEFAULT_EXTREMA_MAX_ITER 200 // default 100
#define _REMEZ_DEFAULT_ROOTS_MAX_ITER 200 // default 100

#define _REMEZ_DEFAULT_RATIONAL_ERR_MAX_CONVGC 200 // default 80
#define _REMEZ_DEFAULT_RATIONAL_ERR_TOLERANCE 0.0001 // default 0.001

#define _REMEZ_DEFAULT_CTRLPT_PERTURBATION 0.005 // default 0.05
#define _REMEZ_DEFAULT_CTRLPT_MAX_PERTURBATION 0.85 // default 0.8

#define _REMEZ_SQRT_DOUBLE_EPSILON 0.00000001490116119384765625 //std::sqrt(std::numeric_limits<double>::epsilon())
#endif

/****
* references:
* boost.org, .\libs\math\include_private\boost\math\tools\remez.hpp
*/

namespace MyMath {
	class Remez {
	private:
		//params
		double (*func)(double);

		bool relative_error; // default true
		bool pinned;

		unsigned int oN;
		unsigned int oD;

		double left_interval;
		double right_interval;

		int skew;

		unsigned int brake;

		unsigned int maxRank;

		unsigned int n;

		//iteration data
		unsigned int iter_count;
		bool sanity;

		double* solution;
		double** matrix;

		RationalFuncContainer* ratFuncC;
		FuncContainer* funcC;

		double* error_func_roots;
		double** error_func_extrema;

		double max_error;
		double max_error_change;
		double max_error_last_change;

		double solution_max_abs_error;
		double solution_max_rel_error;

		bool CheckValidNum(double& x);

		double* GenerateChebyshevKnotsAndBoundaries(); // return double array
		void Initialize();
	public:
		Remez(double (*func)(double)
			, unsigned int oN, unsigned int oD
			, double left_interval, double right_interval
			, int skew);

		Remez(double (*func)(double)
			, bool pinned
			, unsigned int oN, unsigned int oD
			, double left_interval, double right_interval
			, int skew);

		Remez(double (*func)(double)
			, bool relative_error, bool pinned
			, unsigned int oN, unsigned int oD
			, double left_interval, double right_interval
			, int skew);

		~Remez();

		unsigned int GetBrake();
		void SetBrake(unsigned int brake);

		bool Iterate(); // return success

		double GetMaxError();
		double GetCurrentChangeOfMaxError();
		double GetLastChangeOfMaxError();

		double GetSolutionMaxAbsoluteError();
		double GetSolutionMaxRelativeError();

		unsigned int GetIterationCount();
		bool Sanity();

		unsigned int GetErrorFuncRoots(double*& roots);
		unsigned int GetControlPoints(double**& extrema);

		unsigned int GetNumerator(double*& coefficients);
		unsigned int GetDenominator(double*& coefficients);

		double Estimate(double x);
	};
}

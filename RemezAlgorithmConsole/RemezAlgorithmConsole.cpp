// RemezAlgorithmConsole.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//

#include <iostream>
#include <string>
#include <sstream>
#include<conio.h>
#include <Math.h>
#include "Remez.h"

#define PRECISION 52 //23

using namespace std;
using namespace MyMath;

double func(const double x) {
    return exp(x) - 1;
}

int main()
{
    const bool relatErr = true; // relative Error
    const bool pinned = true; // if "pinned" is true, all control points cannot be origin(0,0). The origin(0,0) control point will occur the 0 dividing during LUP step.
    const unsigned int oN = 5; // m, if "pinned" is true, please avoid to use odd rank, they have origin(0,0) Chebyshev Knot
    const unsigned int oD = 5; // n
    const int skew = 0; // in [-100, 100]
    const double a = -1;
    const double b = 1;

    ostringstream out;
    out.precision(PRECISION);

    cout << "Remez Algorithm" << endl;
    cout << endl;
    cout << "Optimize " << (relatErr ? "Relative Error" : "Absolute Error") << endl;
    if (pinned)
        cout << "Pass through the origin (0,0)" << endl;
    cout << "=========================================================================" << endl;
    cout << "Initialization : " << endl;

    Remez remez = { func, relatErr, pinned, oN, oD, a, b, skew };

    double* error_func_roots;
    unsigned int roots_count = remez.GetErrorFuncRoots(error_func_roots);

    double** controlPoints;
    unsigned int controlPoints_count = remez.GetControlPoints(controlPoints);

    double* numerator;
    unsigned int numerator_count = remez.GetNumerator(numerator);

    double* denominator;
    unsigned int denominator_count = remez.GetDenominator(denominator);

    //Chebyshev Knots
    cout << "Initial Points - Chebyshev Knots (Count : " << to_string(roots_count) << ") : " << endl;
    cout << "Skew = " << to_string(skew) << " % " << endl;
    for (unsigned int i = 0; i < roots_count; i++) {
        out.str("");
        out.clear();
        out << std::fixed << error_func_roots[i];
        //cout << "i = " << to_string(i) << " : ";
        cout << "x = ";
        cout << move(out).str() << endl;
    }
    cout << endl;

    //print numerator
    cout << "Numerator : " << endl;
    if (numerator_count > 0) {
        cout << " ";// << endl;
        out.str("");
        out.clear();
        out << std::fixed << numerator[0];
        if (pinned)
            out << "*x^{1}";
        cout << move(out).str();// << endl;
        for (unsigned int i = 1; i < numerator_count; i++) {
            out.str("");
            out.clear();
            out << std::fixed << numerator[i];
            cout << " + " << move(out).str() << "*x^{" << to_string(pinned ? (i + 1) : i) << "}";// << endl;
        }
        cout << " ";// << endl;
    }
    cout << "\nDenominator : " << endl;
    //cout << " / ";
    //print denominator
    if (denominator_count > 0) {
        cout << " 1 + ";// << endl;
        out.str("");
        out.clear();
        out << std::fixed << denominator[0] << "*x^{1}";
        cout << move(out).str();// << endl;
        for (unsigned int i = 1; i < denominator_count; i++) {
            out.str("");
            out.clear();
            out << std::fixed << denominator[i];
            cout << " + " << move(out).str() << "*x^{" << to_string(i + 1) << "}";// << endl;
        }
        cout << " ";// << endl;
    }

    cout << endl << endl;

    out.str("");
    out.clear();
    out << std::fixed << remez.GetMaxError();
    cout << "Max Error = " << move(out).str() << endl;

    cout << endl;

    cout << "Sanity : " << (remez.Sanity() ? "Good" : "Bad") << endl;

    cout << endl;

    cout << "Initialization finished" << endl;
    cout << "=========================================================================" << endl;

    char operKey = 'E';

    //check continue
    cout << "Press Any Key to Iteration; \"C\" for Estimate function; \"E\" for Exit : ";
    operKey = _getch();
    cout << "\n=========================================================================" << endl;

    while (operKey != 'e' && operKey != 'E') {
        if (operKey != 'c' && operKey != 'C') {
            remez.Iterate();

            cout << "Iteration (Count : " << to_string(remez.GetIterationCount()) << ") : " << endl;

            //control points
            cout << "Control Points (Count : " << to_string(controlPoints_count) << ") : " << endl;
            for (unsigned int i = 0; i < controlPoints_count; i++) {
                out.str("");
                out.clear();
                out << std::fixed << controlPoints[i][0];
                //cout << "i = " << to_string(i) << " : ";
                cout << "x = ";
                cout << move(out).str() << endl;
                out.str("");
                out.clear();
                out << std::fixed << controlPoints[i][1];
                //cout << "i = " << to_string(i) << " : ";
                cout << "|f(x)| = ";
                cout << move(out).str() << endl;
            }
            cout << endl;

            //print numerator
            cout << "Numerator : " << endl;
            if (numerator_count > 0) {
                cout << " ";// << endl;
                out.str("");
                out.clear();
                out << std::fixed << numerator[0];
                if (pinned)
                    out << "*x^{1}";
                cout << move(out).str();// << endl;
                for (unsigned int i = 1; i < numerator_count; i++) {
                    out.str("");
                    out.clear();
                    out << std::fixed << numerator[i];
                    cout << " + " << move(out).str() << "*x^{" << to_string(pinned ? (i + 1) : i) << "}";// << endl;
                }
                cout << " ";// << endl;
            }
            cout << "\nDenominator : " << endl;
            //cout << " / ";
            //print denominator
            if (denominator_count > 0) {
                cout << " 1 + ";// << endl;
                out.str("");
                out.clear();
                out << std::fixed << denominator[0] << "*x^{1}";
                cout << move(out).str();// << endl;
                for (unsigned int i = 1; i < denominator_count; i++) {
                    out.str("");
                    out.clear();
                    out << std::fixed << denominator[i];
                    cout << " + " << move(out).str() << "*x^{" << to_string(i + 1) << "}";// << endl;
                }
                cout << " ";// << endl;
            }

            cout << endl << endl;

            double max_err, max_err_change, max_err_last_change;
            out.str("");
            out.clear();
            out << std::fixed << (max_err = remez.GetMaxError());
            cout << "Max Error =   " << (max_err < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl;

            out.str("");
            out.clear();
            out << std::fixed << (max_err_change = remez.GetCurrentChangeOfMaxError());
            cout << "Change =      " << (max_err_change < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl;

            out.str("");
            out.clear();
            out << std::fixed << (max_err_last_change = remez.GetLastChangeOfMaxError());
            cout << "Last Change = " << (max_err_last_change < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl;

            double sol_max_abs_err, sol_max_rel_err;
            out.str("");
            out.clear();
            out << std::fixed << (sol_max_abs_err = remez.GetSolutionMaxAbsoluteError());
            cout << "Solution Max Absolute Error = " << (sol_max_abs_err < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl;

            out.str("");
            out.clear();
            out << std::fixed << (sol_max_rel_err = remez.GetSolutionMaxRelativeError());
            cout << "Solution Max Relative Error = " << (sol_max_rel_err < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl;

            out.str("");
            out.clear();
            out << std::fixed << _REMEZ_SQRT_DOUBLE_EPSILON;
            cout << "Sqrt(EPSILON) = " << move(out).str() << endl;

            cout << endl;

            cout << "Sanity : " << (remez.Sanity() ? "Good" : "Bad") << endl;

            cout << "Iteration finished" << endl;
            cout << "=========================================================================" << endl;
        }
        else {
            cout << "Estimate : " << endl;

            //double x = 0.00020011092601223321;
            double x = 0.0;

            cout << "Please input x in [";
            cout << to_string(a);
            cout << ", ";
            cout << to_string(b);
            cout << "] :" << "\n" << " x = ";
            string inputWord = "";
            char temp;
            while (cin.get(temp) && temp != '\n') {
                if ((temp == ' ') && (static_cast<int>(inputWord.length()) != 0)) {
                    x = stof(inputWord);
                    inputWord = "";
                }
                else if (temp >= 0x30 && temp <= 0x39 || temp == 0x2E || temp == 0x2D) // (0-9 .-)
                    inputWord += temp;
            }
            if (static_cast<int>(inputWord.length()) != 0) {
                x = stof(inputWord);
                inputWord = "";
            }

            cout << endl;
            cout << "-----------------------------------------------------------------------------" << endl;
            //estimate function

            double estimAns = remez.Estimate(x);
            double standardAns = func(x);

            double error = estimAns - standardAns;

            out.str("");
            out.clear();
            out << std::fixed << estimAns;

            cout << "Estimate : " << (estimAns < 0 ? "" : " ") << move(out).str() << endl;

            out.str("");
            out.clear();
            out << std::fixed << standardAns;

            cout << "Standard : " << (standardAns < 0 ? "" : " ") << move(out).str() << endl;

            out.str("");
            out.clear();
            out << std::fixed << error;

            cout << "Error :    " << (error < 0 ? "" : " ") << move(out).str() << endl;

            cout << endl << endl;
        }

        //check continue
        cout << "Press Any Key to Iteration; \"C\" for Estimate function; \"E\" for Exit : ";
        operKey = _getch();
        cout << "\n=========================================================================" << endl;
    }

    system("pause");
}

// 執行程式: Ctrl + F5 或 [偵錯] > [啟動但不偵錯] 功能表
// 偵錯程式: F5 或 [偵錯] > [啟動偵錯] 功能表

// 開始使用的提示: 
//   1. 使用 [方案總管] 視窗，新增/管理檔案
//   2. 使用 [Team Explorer] 視窗，連線到原始檔控制
//   3. 使用 [輸出] 視窗，參閱組建輸出與其他訊息
//   4. 使用 [錯誤清單] 視窗，檢視錯誤
//   5. 前往 [專案] > [新增項目]，建立新的程式碼檔案，或是前往 [專案] > [新增現有項目]，將現有程式碼檔案新增至專案
//   6. 之後要再次開啟此專案時，請前往 [檔案] > [開啟] > [專案]，然後選取 .sln 檔案

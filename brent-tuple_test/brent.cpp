#include "pch.h"
#include "brent.h"

extern double b;

constexpr double r8_epsilon = 2.220446049250313E-016;

double gauss(double x)
{
    double a = 10.0;
    double c = 15.0;
    return -a * exp(-(x - b) * (x - b) / (2 * c * c));
}

std::tuple<int, double> local_min(double a, double b, double t, double f(double x), double& x)
{
	double c;
    double d;
    double e;
    double eps;
    double fu;
    double fv;
    double fw;
    double fx;
    double m;
    double p;
    double q;
    double r;
    double sa;
    double sb;
    double t2;
    double tol;
    double u;
    double v;
    double w;

    bool parabola = false;
    bool gold = false;
    std::vector<double> b_array = {};
    std::vector<double> fb_array = {};
    int step = 0;
    int dir = 1;

    c = 0.5 * (3.0 - sqrt(5.0));

    eps = sqrt(r8_epsilon);

    sa = a;
    sb = b;

    x = sa + c * (b - a);
    b_array.emplace_back(x);

    w = x;
    v = w;
    e = 0.0;

    Sleep(15);
    fx = f(x);
    fb_array.emplace_back(fx);

    fw = fx;
    fv = fw;

    //int count = round(x) / 10;
    //while (count > 0) {
    //    printf(" - ");
    //    count --;
    //}

    printf("vol = %lf  f = %lf\n", x, fx);

    //printf("\nx = %lf, fx = %lf\n", x, fx);
    step ++;

    for (;;)
    {
        m = 0.5 * (sa + sb);
        tol = eps * fabs(x) + t;
        t2 = 2.0 * tol;
        //
        //  Check the stopping criterion.
        //
        if (fabs(x - m) <= t2 - 0.5 * (sb - sa))
        {
            break;
        }
        //
        //  Fit a parabola.
        //
        r = 0.0;
        q = r;
        p = q;

        if (tol < fabs(e))
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (0.0 < q)
            {
                p = -p;
            }
            q = fabs(q);
            r = e;
            e = d;
        }

        if (fabs(p) < fabs(0.5 * q * r) &&
            q * (sa - x) < p &&
            p < q * (sb - x))
        {
            //
            //  Take the parabolic interpolation step.
            //
            d = p / q;
            u = x + d;
            //
            //  F must not be evaluated too close to A or B.
            //
            if ((u - sa) < t2 || (sb - u) < t2)
            {
                if (x < m)
                {
                    d = tol;
                }
                else
                {
                    d = -tol;
                }
            }

//            if (parabola)
//            {
//                printf("b_array.back() = %lf\n", b_array.back());
//                printf("b_array.at(b_array.size() - 2)) = %lf\n", b_array.at(b_array.size() - 2));
//                printf("u = %lf\n", u);
//                printf("u_forward = %lf\n", b_array.back());
//                if (b_array.back() < u)  dir = 1;
//                else  dir = -1;
//
//                return std::make_tuple(step, b_array.back(), dir);
//            }
            parabola = true;
        }
            //
            //  A golden-section step.
            //
        else
        {
            if (x < m)
            {
                e = sb - x;
            }
            else
            {
                e = sa - x;
            }
            d = c * e;
            //gold = true;

//            if (step > 4)
//            {
//                double tmp = fabs(b_array.back() - b_array.at(b_array.size() - 2));
//                if (tmp < 8.0)
//                {
//                    printf("b_array.back() = %lf\n", b_array.back());
//                    printf("b_array.at(b_array.size() - 2)) = %lf\n", b_array.at(b_array.size() - 2));
//                    printf("u = %lf\n", u);
//                    printf("u_forward = %lf\n", b_array.back());
//                    return std::make_tuple(step, x, dir);
//                }
//            }
        }
        //
        //  F must not be evaluated too close to X.
        //
        if (tol <= fabs(d))
        {
            u = x + d;
        }
        else if (0.0 < d)
        {
            u = x + tol;
        }
        else
        {
            u = x - tol;
        }

        b_array.emplace_back(u);

        Sleep(15);
        fu = f(u);
        fb_array.emplace_back(fu);
        //count = round(u) / 10;
        //while (count > 0) {
        //    printf(" - ");
        //    count --;
        //}
//        if (u < x)  {
//            printf("vol = 0 \n");
//            step ++;
//        }
        printf("vol = %lf  f = %lf parabola = %d\n", u, fu, parabola);
        //parabola = false, gold = false;
        //printf("\nu = %lf, fu = %lf\n", u, fu);
        step ++;
        //
        //  Update A, B, V, W, and X.
        //
        if (fu <= fx)
        {
            if (u < x)
            {
                sb = x;
            }
            else
            {
                sa = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if (u < x)
            {
                sa = u;
            }
            else
            {
                sb = u;
            }

            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    //fx = -gauss(x);
    return std::make_tuple(step, -fx);
}

double Hill_climbing(double &x, double &fx)
{
    double f = fx;
    Sleep(2);
    int dir = 0;
    int step = 0;
    x = x + 1.0;

    double fb = -gauss(x); 
    if (fb > f) dir = 1;
    else dir = -1;
	
    double tmp_fb = fb;

    while(fb - tmp_fb >= 0.0) {
		printf("vol = %lf  fb = %lf\n", x, -gauss(x));
        tmp_fb = max(tmp_fb, fb);
        x += dir * 1.0;
    	Sleep(2);
        fb = -gauss(x);
        step++;
        if (x > 150 && x < 0) break;
    }
    return tmp_fb;
}
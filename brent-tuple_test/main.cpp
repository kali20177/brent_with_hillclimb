// brent-tuple_test.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include "brent.h"

double b = 9.0;

int main()
{
    double x = 0.0;
    double fx = 0.0;
    int step = 0;
	
    std::cout << "Brent search:" << std::endl;
    std::tie(step, fx) = local_min(0.0, 80.0, 8.0, gauss, x);
    std::cout << "vol = " << x << "\n";
	std::cout << "fx = " << fx << "\n";
    std::cout << "step = " << step << "\n";

    std::cout << "climb hill:" << std::endl;
	double res = Hill_climbing(x, fx);
    std::cout <<  "x = " << x << "\n";
    std::cout <<  "maxfb = " << res << "\n";
    return 0;
}


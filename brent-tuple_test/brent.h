#pragma once

double gauss(double x);

std::tuple<int, double> local_min(double a, double b, double t, double f(double x), double& x);

double Hill_climbing(double &x, double &fx);
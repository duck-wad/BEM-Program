#pragma once
#include <vector>

double Potential(double r, double k, int Cdim);
double Flow(double r, std::vector<double>& dxr, std::vector<double>& normal, int Cdim);
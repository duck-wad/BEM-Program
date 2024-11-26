#pragma once

#include <vector>

const double PI = 3.14159265359;

const double LOW_TOL = 1e-8;
const double HIGH_TOL = 1e-12;

std::vector<std::vector<double>>& operator*=(std::vector<std::vector<double>>& matrix, double constant);

void PrintMatrix(const std::vector<std::vector<double>>& matrix);

int Kronecker(int i, int j);
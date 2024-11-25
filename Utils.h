#pragma once

#include <vector>

const double PI = 3.14159265359;

std::vector<std::vector<double>>& operator*=(std::vector<std::vector<double>>& matrix, double constant);

void PrintMatrix(const std::vector<std::vector<double>>& matrix);
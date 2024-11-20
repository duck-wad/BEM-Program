#pragma once

#include <vector>
#include "Enum.h"

double DotProduct(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> CrossProduct(const std::vector<double>& a, const std::vector<double>& b);
double VectorLength(const std::vector<double>& a);

void NormalizeVector(std::vector<double>& a, const double L);
void NormalJacobian(std::vector<double>& v3, double& Jac, double xsi, double eta, ElementType type, 
	int nodes, const std::vector<int>& inci, const std::vector<std::vector<double>>& coords);
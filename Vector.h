#pragma once

#include <vector>

double DotProduct(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> CrossProduct(const std::vector<double>& a, const std::vector<double>& b);
double VectorLength(const std::vector<double>& a);

std::vector<double> NormalizeVector(const std::vector<double>& a, const double L);
void NormalJacobian(std::vector<double>& v3, double& Jac, double xsi, double eta, int ldim, 
	int nodes, const std::vector<int>& inci, const std::vector<std::vector<double>>& coords);
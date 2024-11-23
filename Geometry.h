#pragma once

#include <vector>

void SerendipityFunction(std::vector<double>& Ni, int ldim,
	double xsi, double eta, int nodes, const std::vector<int>& inci);

void SerendipityDerivative(std::vector<std::vector<double>>& DNi, int ldim,
	double xsi, double eta, int nodes, const std::vector<int>& inci);

void ComputeArea();
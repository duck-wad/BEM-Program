#pragma once

#include <vector>

#include "Enum.h"

void SerendipityFunction(std::vector<double>& Ni, ElementType type,
	double xsi, double eta, int nodes, const std::vector<int>& inci);

void SerendipityDerivative(std::vector<std::vector<double>>& DNi, ElementType type,
	double xsi, double eta, int nodes, const std::vector<int>& inci);
#pragma once

#include <vector>

//receives an integration order, outputs vector of two vectors
//first vector contains the gauss point coordinates, second contains the weights
void GaussPoints(int order, std::vector<double>& gaussPoints, std::vector<double>& weights);
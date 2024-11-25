#include <iostream>
#include <cmath>

#include "Vector.h"
#include "Utils.h"

//receives radial distance and permeability/conductivity, outputs the potential
double Potential(double r, double k, int Cdim) {
	double U = 0.0;

	if (Cdim == 2){
		U = 1 / (2 * PI * k) * log(1 / r);
	}
	else if (Cdim == 3) {
		U = 1 / (4 * PI * r * k);
	}
	else {
		throw std::invalid_argument("Dimension must be 2 or 3");
	}
	return U;
}

//receives radial distance (magnitude), distance vector, and normal vector to boundary, outputs the flow through the boundary
double Flow(double r, std::vector<double>& dxr, std::vector<double>& normal, int Cdim) {
	double T = 0.0;

	//normal and dxr are both expected to be unit vectors prior to input to this function
	if (Cdim == 2) {
		T = DotProduct(dxr, normal) / (2 * PI * r);
	}
	else if (Cdim == 3) {
		T = DotProduct(dxr, normal) / (4 * PI * r * r);
	}
	else {
		throw std::invalid_argument("Dimension must be 2 or 3");
	}
	return T;
}
#include <iostream>
#include <cmath>
#include <cassert>

#include "Kernels.h"
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

//receives radial distance (magnitude), distance unit vector, and normal vector to boundary, outputs the flow through the boundary
double Flow(double r, const std::vector<double>& dxr, const std::vector<double>& normal, int Cdim) {
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

//receives radial distance, distance unit vector, Young's modulus, Poisson ratio. outputs rank 2 or 3 matrix (depending on dimension) for the displacements on boundary at point Q due to unit point load at P in the domain
void KelvinDisplacement(std::vector<std::vector<double>>& UK, const std::vector<double>& dxr, double r, double E, double nu, int Cdim) {

	double G = E / (2.0 * (1.0 + nu));
	double c;
	double c1 = 3.0 - 4.0 * nu;

	if (Cdim == 2) {
		c = 1.0 / (8.0 * PI * G * (1.0 - nu));
		UK.resize(2, std::vector<double>(2, 0.0));
		for (size_t i = 0; i < Cdim; i++) {
			for (size_t j = 0; j < Cdim; j++) {
				UK[i][j] = c * (c1 * log(1 / r) * Kronecker(static_cast<int>(i), static_cast<int>(j)) + dxr[i] * dxr[j]);
			}
		}
	}
	else if (Cdim == 3) {
		c = 1.0 / (16.0 * PI * G * (1 - nu));
		UK.resize(3, std::vector<double>(3, 0.0));
		for (size_t i = 0; i < Cdim; i++) {
			for (size_t j = 0; j < Cdim; j++) {
				UK[i][j] = c * (c1 * Kronecker(static_cast<int>(i), static_cast<int>(j)) + dxr[i] * dxr[j]);
			}
		}
	}
	else
		throw std::invalid_argument("Dimension must be 2 or 3");
}

//output rank 2 or 3 matrix for traction on boundary at Q due to unit point load at P. note that the dxr and normal vectors are expected to be normalized unit vectors
void KelvinTraction(std::vector<std::vector<double>>& TK, const std::vector<double>& dxr, const std::vector<double>& normal, double r, double nu, int Cdim) {

	double c2;
	double c3 = 1.0 - 2.0 * nu;
	int temp;

	//assume normal and dxr are unit length
	assert(VectorLength(normal) - 1.0 <= LOW_TOL);
	assert(VectorLength(dxr) - 1.0 <= LOW_TOL);

	double costheta = DotProduct(normal, dxr);
	
	if (Cdim == 2) {
		c2 = 1.0 / (4.0 * PI * (1.0 - nu));
		temp = 2;
		TK.resize(2, std::vector<double>(2, 0.0));
	}
	else if (Cdim == 3) {
		c2 = 1.0 / (8.0 * PI * (1.0 - nu));
		temp = 3;
		TK.resize(3, std::vector<double>(3, 0.0));
	}
	else
		throw std::invalid_argument("Dimension must be 2 or 3");

	for (size_t i = 0; i < Cdim; i++) {
		for (size_t j = 0; j < Cdim; j++) {
			TK[i][j] = c2 / (r * r) * (c3 * Kronecker(static_cast<int>(i), static_cast<int>(j)) + temp * dxr[i] * dxr[j] * costheta - c3 * (1 - Kronecker(static_cast<int>(i), static_cast<int>(j))) * (normal[j] * dxr[i] - normal[i] * dxr[j]));
		}
	}
}
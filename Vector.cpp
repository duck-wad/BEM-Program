#include <iostream>b
#include <cassert>

#include "Vector.h"
#include "Serendipity.h";

double DotProduct(const std::vector<double>& a, const std::vector<double>& b) {

	assert(a.size() == b.size() && "Vectors must be same size to dot product");

	double result = 0.;

	for (size_t i = 0; i < a.size(); i++) {
		result += a[i] * b[i];
	}
	return result;
}

std::vector<double> CrossProduct(const std::vector<double>& a, const std::vector<double>& b) {

	assert(a.size() && b.size() == 3 && "Vectors are not the right dimension for cross product");

	std::vector<double> result(3, 0);

	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = -1*(a[0] * b[2] - a[2] * b[0]);
	result[2] = a[0] * b[1] - a[1] * b[0];
	return result;
}

double VectorLength(const std::vector<double>& a) {

	double L = 0.;

	for (size_t i = 0; i < a.size(); i++) {
		L += a[i] * a[i];
	}

	return sqrt(L);
}

//normalize vector by its length
void NormalizeVector(std::vector<double>& a, const double L) {

	assert(L != 0 && "Cannot normalize vector with zero length");

	for (size_t i = 0; i < a.size(); i++) {
		a[i] /= L;
	}
}

//input point coordinates, output the Jacobian and normal vector to the point
void NormalJacobian(std::vector<double>& v3, double& Jac, double xsi, double eta, ElementType type,
	int nodes, const std::vector<int>& inci, const std::vector<std::vector<double>>& coords) {

	int Cdim = (type == QUAD) ? 3 : 2;

	std::vector<std::vector<double>> DNi(2, std::vector<double>(nodes, 0.0));
	//initialize tangent vectors in xsi (v1) and eta (v2) direction
	std::vector<double> v1(Cdim, 0.), v2(Cdim, 0.);

	SerendipityDerivative(DNi, type, xsi, eta, nodes, inci);

	//compute the tangent vectors in xsi and eta direction
	for (size_t i = 0; i < Cdim; i++) {
		v1[i] = DotProduct(DNi[0], coords[i]);
		if (type == QUAD) {
			v2[i] = DotProduct(DNi[1], coords[i]);
		}
	}
	//compute the normal vector 
	if (type == LINE) {
		v3[0] = v1[1];
		v3[1] = -1 * v1[0];
	}
	else if (type == QUAD) {
		v3 = CrossProduct(v1, v2);
	}

	Jac = VectorLength(v3);

	NormalizeVector(v3, Jac);

	//remove negative signs from zero values
	for (size_t i = 0; i < v3.size(); i++) {
		if (v3[i] <= DBL_EPSILON) {
			v3[i] = std::fabs(v3[i]);
		}
	}
}
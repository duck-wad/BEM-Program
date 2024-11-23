#include <iostream>
#include <cassert>

#include "Integration.h"

void GaussPoints(int order, std::vector<double>& gaussPoints, std::vector<double>& weights) {

	gaussPoints.resize(order);
	weights.resize(order);

	switch (order) {
	case 1:
		gaussPoints[0] = 0.0;

		weights[0] = 2.0;
		break;
	case 2:
		gaussPoints[0] = 0.577350269;
		gaussPoints[1] = -0.577350269;

		weights[0] = weights[1] = 1.0;
		break;
	case 3:
		gaussPoints[0] = 0.774596669;
		gaussPoints[1] = 0.0;
		gaussPoints[2] = -0.774596669;

		weights[0] = weights[2] = 0.555555555;
		weights[1] = 0.888888888;
		break;
	case 4:
		gaussPoints[0] = 0.861136311;
		gaussPoints[1] = 0.339981043;
		gaussPoints[2] = -0.339981043;
		gaussPoints[3] = -0.861136311;

		weights[0] = weights[3] = 0.347854845;
		weights[1] = weights[2] = 0.652145154;
		break;
	case 5:
		gaussPoints[0] = 0.906179845;
		gaussPoints[1] = 0.538469310;
		gaussPoints[2] = 0.0;
		gaussPoints[3] = -0.538469310;
		gaussPoints[4] = -0.906179845;

		weights[0] = weights[4] = 0.236926885;
		weights[1] = weights[3] = 0.478628670;
		weights[2] = 0.568888888;
		break;
	case 6:
		gaussPoints[0] = 0.932469514;
		gaussPoints[1] = 0.661209386;
		gaussPoints[2] = 0.238619186;
		gaussPoints[3] = -0.238619186;
		gaussPoints[4] = -0.661209386;
		gaussPoints[5] = -0.932469514;

		weights[0] = weights[5] = 0.171324492;
		weights[1] = weights[4] = 0.360761573;
		weights[2] = weights[3] = 0.467913934;
		break;
	case 7:
		gaussPoints[0] = 0.949107912;
		gaussPoints[1] = 0.741531185;
		gaussPoints[2] = 0.405845151;
		gaussPoints[3] = 0.0;
		gaussPoints[4] = -0.405845151;
		gaussPoints[5] = -0.741531185;
		gaussPoints[6] = -0.949107912;

		weights[0] = 0.129484966;
		weights[1] = 0.279705391;
		weights[2] = 0.381830050;
		weights[3] = 0.417959183;
		weights[4] = weights[2];
		weights[5] = weights[1];
		weights[6] = weights[0];
		break;
	case 8:
		gaussPoints[0] = 0.960289856;
		gaussPoints[1] = 0.796666477;
		gaussPoints[2] = 0.525532409;
		gaussPoints[3] = 0.183434642;
		gaussPoints[4] = -0.183434642;
		gaussPoints[5] = -0.525532409;
		gaussPoints[6] = -0.796666477;
		gaussPoints[7] = -0.960289856;

		weights[0] = 0.101228536;
		weights[1] = 0.222381034;
		weights[2] = 0.313706645;
		weights[3] = 0.362683783;
		weights[4] = weights[3];
		weights[5] = weights[2];
		weights[6] = weights[1];
		weights[7] = weights[0];
		break;
	default:
		throw std::invalid_argument("Invalid order for gauss quadrature");
	}
}
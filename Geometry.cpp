#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>

#include "Geometry.h"

//function receives a vector with shape function values, xsi and eta local coordinates, number of dimensions, number of nodes, and vector representing node indices
//return the value of the shape functions at the given local coordinates
//incidence vector contains the node numbers, if one of them is zero then the node is not active (missing)
void SerendipityFunction(std::vector<double>& Ni, int ldim,
	double xsi, double eta, int nodes, const std::vector<int>& inci) {

	assert(static_cast<int>(inci.size() == nodes) && "Error: incidence vector must size must match number of nodes");

	//temporary variables to simplify calculations
	double mxs, pxs, met, pet;

	if (ldim == 1) {

		assert((nodes == 2 || nodes == 3) && "Error: 1D elements must have 2 or 3 nodes");

		Ni[0] = 0.5 * (1. - xsi);
		Ni[1] = 0.5 * (1. + xsi);

		if (nodes == 2) return;

		Ni[2] = 1. - xsi * xsi;
		Ni[0] = Ni[0] - 0.5 * Ni[2];
		Ni[1] = Ni[1] - 0.5 * Ni[2];
	}
	else if (ldim == 2) {

		assert((nodes == 4 || nodes == 8) && "Error: quadrilateral elements must have 4 or 8 nodes");

		mxs = 1.0 - xsi;
		pxs = 1.0 + xsi;
		met = 1.0 - eta;
		pet = 1.0 + eta;

		Ni[0] = 0.5 * mxs * 0.5 * met;
		Ni[1] = 0.5 * pxs * 0.5 * met;
		Ni[2] = 0.5 * pxs * 0.5 * pet;
		Ni[3] = 0.5 * mxs * 0.5 * pet;

		if (nodes == 4) return;

		if (inci[4] > 0) {
			Ni[4] = 0.5 * (1. - xsi * xsi) * met;
			Ni[0] = Ni[0] - 0.5 * Ni[4];
			Ni[1] = Ni[1] - 0.5 * Ni[4];
		}
		if (inci[5] > 0) {
			Ni[5] = 0.5 * (1. - eta * eta) * pxs;
			Ni[1] = Ni[1] - 0.5 * Ni[5];
			Ni[2] = Ni[2] - 0.5 * Ni[5];
		}
		if (inci[6] > 0) {
			Ni[6] = 0.5 * (1. - xsi * xsi) * pet;
			Ni[2] = Ni[2] - 0.5 * Ni[6];
			Ni[3] = Ni[3] - 0.5 * Ni[6];
		}
		if (inci[7] > 0) {
			Ni[7] = 0.5 * (1. - eta * eta) * mxs;
			Ni[0] = Ni[0] - 0.5 * Ni[7];
			Ni[3] = Ni[3] - 0.5 * Ni[7];
		}
	}
	else {
		std::cout << "Error: not a valid element type" << std::endl;
	}
}

//function receives vector of vectors for the derivatives
//if element is 1D then the second column of DNi is filled with zeros
//if element is 2D then the first column contains partial derivative wrt 
void SerendipityDerivative(std::vector<std::vector<double>>& DNi, int ldim,
	double xsi, double eta, int nodes, const std::vector<int>& inci) {

	assert(static_cast<int>(inci.size() == nodes) && "Error: incidence vector must size must match number of nodes");

	DNi.resize(2, std::vector<double>(nodes, 0.0));

	//temporary variables to simplify calculations
	double mxs, pxs, met, pet;

	if (ldim == 1) {
		assert((nodes == 2 || nodes == 3) && "Error: 1D elements must have 2 or 3 nodes");

		DNi[0][0] = -0.5;
		DNi[0][1] = 0.5;

		if (nodes == 3) {
			DNi[0][2] = -2. * xsi;
			DNi[0][0] = DNi[0][0] - 0.5 * DNi[0][2];
			DNi[0][1] = DNi[0][1] - 0.5 * DNi[0][2];

			DNi[1][0] = DNi[1][1] = DNi[1][2] = 0.0; //derivatives wrt eta are zero
		}
		else {
			DNi[1][0] = DNi[1][1] = 0.0;
		}
	}
	else if (ldim == 2) {
		//partial derivatives for 2D in xsi and eta direction
		assert((nodes == 4 || nodes == 8) && "Error: quadrilateral elements must have 4 or 8 nodes");

		mxs = 1. - xsi;
		pxs = 1. + xsi;
		met = 1. - eta;
		pet = 1. + eta;

		DNi[0][0] = -0.25 * met;
		DNi[1][0] = -0.25 * mxs;
		DNi[0][1] = 0.25 * met;
		DNi[1][1] = -0.25 * pxs;
		DNi[0][2] = 0.25 * pet;
		DNi[1][2] = 0.25 * pxs;
		DNi[0][3] = -0.25 * pet;
		DNi[1][3] = 0.25 * mxs;

		if (nodes == 4) return;
		if (inci[4] > 0) {
			DNi[0][4] = -xsi * met;
			DNi[1][4] = -0.5 * (1.0 - xsi * xsi);
			//update corner nodes
			DNi[0][0] -= 0.5 * DNi[0][4];
			DNi[1][0] -= 0.5 * DNi[1][4];
			DNi[0][1] -= 0.5 * DNi[0][4];
			DNi[1][1] -= 0.5 * DNi[1][4];
		}
		if (inci[5] > 0) {
			DNi[0][5] = 0.5 * (1. - eta * eta);
			DNi[1][5] = -eta * pxs;

			DNi[0][1] -= 0.5 * DNi[0][5];
			DNi[1][1] -= 0.5 * DNi[1][5];
			DNi[0][2] -= 0.5 * DNi[0][5];
			DNi[1][2] -= 0.5 * DNi[1][5];
		}
		if (inci[6] > 0) {
			DNi[0][6] = -xsi * pet;
			DNi[1][6] = 0.5 * (1. - xsi * xsi);

			DNi[0][2] -= 0.5 * DNi[0][6];
			DNi[1][2] -= 0.5 * DNi[1][6];
			DNi[0][3] -= 0.5 * DNi[0][6];
			DNi[1][3] -= 0.5 * DNi[1][6];
		}
		if (inci[7] > 0) {
			DNi[0][7] = -0.5 * (1. - eta * eta);
			DNi[1][7] = -eta * mxs;

			DNi[0][0] -= 0.5 * DNi[0][7];
			DNi[1][0] -= 0.5 * DNi[1][7];
			DNi[0][3] -= 0.5 * DNi[0][7];
			DNi[1][3] -= 0.5 * DNi[1][7];
		}
	}
	else {
		std::cout << "Error: not a valid element type" << std::endl;
	}
}

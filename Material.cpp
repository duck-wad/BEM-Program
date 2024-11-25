#include <iostream>

#include "Material.h"
#include "Utils.h"

//computes the isotropic D-matrix to relate the strain vector to stress vector
void DMatrix(double E, double nu, int Cdim, std::vector<std::vector<double>>& D) {
	double c1 = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
	double c2 = nu / (1.0 - nu);
	double G = E / (2.0 * (1 + nu));

	if (Cdim == 2) {
		D.resize(3, std::vector<double>(3, 0.0));
		D[0][0] = 1.0;
		D[1][1] = 1.0;
		D[2][2] = G / c1;
		D[1][0] = c2;
		D[0][1] = c2;
	}
	else if (Cdim == 3) {
		D.resize(6, std::vector<double>(6, 0.0));
		D[0][0] = 1.0;
		D[1][1] = 1.0;
		D[2][2] = 1.0;
		D[1][0] = c2;
		D[0][1] = c2;
		D[0][2] = c2;
		D[2][0] = c2;
		D[1][2] = c2;
		D[2][1] = c2;
		D[3][3] = G / c1;
		D[4][4] = G / c1;
		D[5][5] = G / c1;
	}
	D *= c1;
}
#include <iostream>
#include <vector>
#include <cassert>

//#include "Enum.h"
#include "Serendipity.h"
#include "Vector.h"


int main() {
	
	ElementType type = QUAD;
	int nodes = 8;
	double xsi = 0.1, eta = 0.2;
	std::vector<int> inci = { 1,2,3,4,5,6,7,8 };
	std::vector<double> Ni(nodes, 0.0);
	std::vector<std::vector<double >> DNi;

	SerendipityFunction(Ni, type, xsi, eta, nodes, inci);
	SerendipityDerivative(DNi, type, xsi, eta, nodes, inci);

	for (int i = 0; i < DNi.size(); ++i) {
		std::cout << "Node " << i + 1 << " derivatives:\n";
		for (int j = 0; j < DNi[i].size(); ++j) {
			if (j == 0) {
				std::cout << "  dNi/dxsi = " << DNi[i][j] << "\n";
			}
			else if (j == 1) {
				std::cout << "  dNi/deta = " << DNi[i][j] << "\n";
			}
		}
		std::cout << "\n";
	}

	/*
	std::vector<std::vector<double>> coords = {
		{0.0, 1.0},  // x-coordinates of the nodes
		{0.0, 0.0}   // y-coordinates of the nodes
	};
	std::vector<int> inci = { 1, 2 };
	ElementType type = LINE;
	int nodes = 2;
	double xsi = 0.0, eta = 0.0;
	std::vector<double> v3(2, 0.0);
	double Jac = 0.0;

	NormalJacobian(v3, Jac, xsi, eta, type, nodes, inci, coords);

	// Expected values:
	// Normal vector (v3): {0.0, 1.0}
	// Jacobian (Jac): 1.0
	std::cout << "LINE Test 1: Normal = {" << v3[0] << ", " << v3[1] << "}, Jacobian = " << Jac << std::endl;
	*/
	
	return 0;
}
#include <iostream>
#include <vector>

#include "Utils.h"

#include "TestCases.h"
#include "Vector.h"
#include "Serendipity.h"

void TestNormalJacobian(int caseNum) {

	std::vector<double> v3;
	double Jac = 0.0;
	double xsi = 0.0; // Center of the curve
	double eta = 0.0; // Not used for LINE elements
	ElementType type;
	int nodes = 2;
	std::vector<int> inci;
	std::vector<std::vector<double>> coords;

	switch (caseNum) {

	case 1:
		//test 2-node line element flat
	case 2:
		//test 2-node line element vertical
	case 3: 
		//test 2-node line element diagonal
	case 4: 
		//test 3-node line element flat
	case 5:
		//test 3-node line element curved like parabola about y-axis
		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0; // Center of the curve
		eta = 0.0; // Not used for LINE elements
		type = LINE;
		nodes = 3;
		inci.resize(3);
		inci = { 1, 2, 3 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 1.0, 0.5 };	// x-coordinates
		coords[1] = { 0.0, 0.0, 1.0 };	// y-coordinates

		NormalJacobian(v3, Jac, xsi, eta, type, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 6: 
		//test 3-node line element curved like parabola about the x-axis
	case 7:
		//4-node quad element flat in xy plane. expecting a vertical normal vector
		v3.resize(3, 0.0);
		Jac = 0.0;
		xsi = 0.0;
		eta = 0.0;
		type = QUAD;
		nodes = 4;
		inci.resize(4);
		inci = { 1,2,3,4 };
		coords.resize(3, std::vector<double>(nodes));
		coords[0] = { 0.0, 1.0, 1.0, 0.0 };
		coords[1] = { 0.0, 0.0, 1.0, 1.0 };
		coords[2] = { 0.0, 0.0, 0.0, 0.0 };

		NormalJacobian(v3, Jac, xsi, eta, type, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << ", " << v3[2] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;
	default:
		std::cout << "Invalid case number" << std::endl;
		return;
	}
}

void TestSerendipity(int caseNum) {

	ElementType type;
	int nodes = 2;
	double xsi = 0.0;
	double eta = 0.0;
	std::vector<int> inci;
	std::vector<double> Ni;
	std::vector<std::vector<double>> DNi;

	switch (caseNum) {
	case 1:
		//2-node line element
	case 2: 
		//3-node line element
	case 3:
		//4-node quad element
	case 4:
		//8-node quad element
		type = QUAD;
		nodes = 8;
		xsi = 0.1;
		eta = 0.2;
		inci = { 1,2,3,4,5,6,7,8 };
		Ni.resize(8, 0.0);
		DNi.resize(2, std::vector<double>(8,0.0));

		SerendipityFunction(Ni, type, xsi, eta, nodes, inci);
		SerendipityDerivative(DNi, type, xsi, eta, nodes, inci);
		
		std::cout << "Element shape functions:" << std::endl;
		for (int i = 0; i < nodes; ++i) {
			std::cout << "  Ni[" << i + 1 << "] = " << Ni[i] << std::endl;
		}
		std::cout << std::endl;

		for (size_t i = 0; i < DNi[0].size(); ++i) {
			std::cout << "Node " << i + 1 << " derivatives:" << std::endl;
			std::cout << "  dNi/dxsi = " << DNi[0][i] << std::endl;
			std::cout << "  dNi/deta = " << DNi[1][i] << std::endl;
		}
		std::cout << std::endl;
		break;
	}

}


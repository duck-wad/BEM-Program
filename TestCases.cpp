#include <iostream>
#include <vector>

#include "TestCases.h"
#include "Vector.h"
#include "Geometry.h"
#include "Integration.h"
#include "Material.h"
#include "Utils.h"

void TestNormalJacobian(int caseNum) {

	std::vector<double> v3;
	double Jac = 0.0;
	double xsi = 0.0; // Center of the curve
	double eta = 0.0; // Not used for LINE elements
	int ldim;
	int nodes = 2;
	std::vector<int> inci;
	std::vector<std::vector<double>> coords;

	switch (caseNum) {

	case 1:
		//test 2-node line element flat
		v3.resize(2, 0.0);
		Jac = 0.0;          
		xsi = 0.0;        
		eta = 0.0;          
		ldim = 1;  
		nodes = 2;
		inci.resize(2);
		inci = { 1,2 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 5.0 };
		coords[1] = { 0.0, 0.0 };

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 2:
		//test 2-node line element vertical

		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0;
		eta = 0.0;
		ldim = 1;
		nodes = 2;
		inci.resize(2);
		inci = { 1,2 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 0.0 };
		coords[1] = { 0.0, 1.0 };

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 3: 
		//test 2-node line element diagonal

		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0;
		eta = 0.0;
		ldim = 1;
		nodes = 2;
		inci.resize(2);
		inci = { 1,2 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 5.0 };
		coords[1] = { 0.0, 5.0 };

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 4: 
		//test 3-node line element flat

		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0; // Center of the curve
		eta = 0.0; // Not used for LINE elements
		ldim = 1;
		nodes = 3;
		inci.resize(3);
		inci = { 1, 2, 3 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 1.0, 0.5 };	// x-coordinates
		coords[1] = { 0.0, 1.0, 0.5};	// y-coordinates

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 5:
		//test 3-node line element curved like parabola about y-axis
		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0; // Center of the curve
		eta = 0.0; // Not used for LINE elements
		ldim = 1;
		nodes = 3;
		inci.resize(3);
		inci = { 1, 2, 3 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 1.0, 0.5 };	// x-coordinates
		coords[1] = { 0.0, 0.0, 1.0 };	// y-coordinates

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 6: 
		//test 3-node line element curved like parabola about the x-axis
		v3.resize(2, 0.0);
		Jac = 0.0;
		xsi = 0.0; // Center of the curve
		eta = 0.0; // Not used for LINE elements
		ldim = 1;
		nodes = 3;
		inci.resize(3);
		inci = { 1, 2, 3 };
		coords.resize(2, std::vector<double>(nodes));
		coords[0] = { 0.0, 0.0, 1.0 };	// x-coordinates
		coords[1] = { 0.0, 1.0, 0.5 };	// y-coordinates

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;

		break;
	case 7:
		//4-node quad element flat in xy plane. expecting a vertical normal vector
		v3.resize(3, 0.0);
		Jac = 0.0;
		xsi = 0.0;
		eta = 0.0;
		ldim = 2;
		nodes = 4;
		inci.resize(4);
		inci = { 1,2,3,4 };
		coords.resize(3, std::vector<double>(nodes));
		coords[0] = { 0.0, 1.0, 1.0, 0.0 };
		coords[1] = { 0.0, 0.0, 1.0, 1.0 };
		coords[2] = { 0.0, 0.0, 0.0, 0.0 };

		NormalJacobian(v3, Jac, xsi, eta, ldim, nodes, inci, coords);

		std::cout << "Test Case: " << caseNum << std::endl;;
		std::cout << "Normal Vector: [" << v3[0] << ", " << v3[1] << ", " << v3[2] << "]" << std::endl;
		std::cout << "Jacobian: " << Jac << std::endl;
	default:
		std::cout << "Invalid case number" << std::endl;
		return;
	}
}

void TestGaussPoints(int order) {
	try {
		// Call your GaussPoints function to get the points and weights
		std::vector<double> gaussPoints;
		std::vector<double> weights;
		GaussPoints(order, gaussPoints, weights);

		std::cout << "Gauss Points and Weights for Order " << order << ":\n";
		std::cout << "-----------------------------------\n";
		std::cout << "Point\t\tWeight\n";

		for (int i = 0; i < order; ++i) {
			std::cout << gaussPoints[i] << "\t\t" << weights[i] << "\n";
		}
		std::cout << "-----------------------------------\n";
	}
	catch (const std::invalid_argument& e) {
		std::cerr << "Error: " << e.what() << "\n";
	}
}

void TestDMatrix() {
	double E = 200000;
	double nu = 0.3;
	int Cdim = 2;
	std::vector<std::vector<double>> mat;
	DMatrix(E, nu, Cdim, mat);
	PrintMatrix(mat);
}
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

#include "Compute.h"
#include "Integration.h"
#include "Vector.h"

void ComputeArea() {
	//ldim is dimension of the boundary element (1D, 2D)
	//lnodes is the number of nodes per element
	//noelem is number of elements
	//intord is the integration order (for Gaussian quadrature)
	//maxnod is the largest node number
	//Cdim is the dimension of analysis (2D for line, 3D for quad)
	int ldim, lnodes, noelem, intord, maxnod, Cdim = 0;
	//vector containing all element indices (global incidence array)
	std::vector<std::vector<int>> inciG;
	//vector containing a single element's indices
	std::vector<int> inci;
	//vector containing coordinates of all element nodes
	std::vector<std::vector<double>> corG;
	//vector containing coordinates of a single element's nodes
	std::vector<std::vector<double>> cor;
	//vector containing normal to the element
	std::vector<double> v3;
	//vector containing gauss points for the integration order
	std::vector<double> Gcor;
	//vector containing corresponding gauss point weights
	std::vector<double> Wi;
	double Jac = 0.0, xsi = 0.0, eta = 0.0, Area = 0.0;

	std::string junk;

	std::ifstream infile("INPUT.DAT");
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
		return;
	}

	infile >> junk >> ldim >> junk >> lnodes >> junk >> noelem >> junk >> intord;
	Cdim = ldim + 1;

	std::cout << "Element dimension = " << ldim << "\n";
	std::cout << "No. of element nodes = " << lnodes << "\n";
	std::cout << "Number of elements = " << noelem << "\n";
	std::cout << "Integration order = " << intord << "\n";

	v3.resize(Cdim);
	inciG.resize(noelem, std::vector<int>(lnodes));

	//reading the global indices
	maxnod = 0;
	for (size_t i = 0; i < noelem; i++) {
		infile >> junk >> junk >> junk;
		for (size_t j = 0; j < lnodes; j++) {
			infile >> inciG[i][j];
			if (inciG[i][j] > maxnod) {
				maxnod = inciG[i][j];
			}
		}
	}

	//global coordinate container has Cdim vectors, one for each direction
	//size is max node number + 1 to have a "node 0" which is set to 0.0
	//this is to handle any cases where the node is missing (i.e node number set to 0)
	//indexing will then start from 1 instead of 0 and go up to max node
	corG.resize(Cdim, std::vector<double>(maxnod + 1, 0.0));

	//initialize node 0 to 0.0 value
	for (size_t i = 0; i < Cdim; i++) {
		corG[i][0] = 0.0;
	}

	//reading the coordinates for each node
	for (size_t i = 1; i <= maxnod; i++) {
		infile >> junk >> junk >> junk;
		for (size_t j = 0; j < Cdim; j++) {
			infile >> corG[j][i];
		}
	}

	inci.resize(lnodes);
	cor.resize(Cdim, std::vector<double>(lnodes));
	Gcor.resize(intord);
	Wi.resize(intord);

	GaussPoints(intord, Gcor, Wi);
	
	//loop over all the elements and add each element area/length to the total
	for (size_t elem = 0; elem < noelem; elem++) {
		inci = inciG[elem];
		for (size_t i = 0; i < Cdim; i++) {
			for (size_t j = 0; j < lnodes; j++) {
				cor[i][j] = corG[i][inci[j]];
			}
		}
		if (ldim == 1) {
			//for 1 dimensional problem, "area" is the length of the boundary
			for (size_t i = 0; i < intord; i++) {
				xsi = Gcor[i];
				NormalJacobian(v3, Jac, xsi, eta, ldim, lnodes, inci, cor);
				Area += Jac * Wi[i];
			}
		}
		else if (ldim == 2) {
			for (size_t i = 0; i < intord; i++) {
				for (size_t j = 0; j < intord; j++) {
					//for 2D elements, need to sum J times Wi for every combination of xsi and eta
					xsi = Gcor[i];
					eta = Gcor[j];
					NormalJacobian(v3, Jac, xsi, eta, ldim, lnodes, inci, cor);
					Area += Jac * Wi[i] * Wi[j];
				}
			}
		}
		else
			std::cout << "Invalid dimension" << std::endl;
	}

	if (ldim == 1) {
		std::cout << "Length = " << Area << std::endl;
	}
	else {
		std::cout << "Area = " << Area << std::endl;
	}
}
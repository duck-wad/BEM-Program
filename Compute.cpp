#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Compute.h"
#include "Integration.h"
#include "Vector.h"
#include "Utils.h"

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

	std::ifstream infile("input//INPUT_AREA.txt");
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
	//area per element is the jacobian at each 
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

//computes the heat flow past an isolated cylinder using the direct BE method
//input file gives flow, conductivity, radius of cylinder and number of segments
//function discretizes the boundary into nseg, 
void FlowAroundCylinder() {

	std::cout << "Heat flow past cylinder using the direct BE method." << std::endl;

	double q, k, radius, radiusOuter, deltaTheta;
	int numseg;
	//kernel matrices
	std::vector<std::vector<double>> deltaT;
	std::vector<std::vector<double>> deltaU;
	//vector with segment center potential
	std::vector<double> u;
	//vector for deltaU times t0 (applied flow)
	std::vector<double> F;
	//applied flow i.e boundary condition
	std::vector<double> t0;
	//start coordA and end coordB coordinates of each segment
	//first row contains the x coord, second row contains the y coord
	//each column correspond with a point
	std::vector<std::vector<double>> coordA, coordB;
	//coordinate of the center point Pi
	std::vector<std::vector<double>> coordPi;
	//temp tangent vector for each segment
	std::vector<double> Vt(2);
	//temp normal unit vector for each segment
	std::vector<double> Vn(2);
	//temp vectors pointing from Pi to A and B of segment
	std::vector<double> vrA(2);
	std::vector<double> vrB(2);
	//segment length
	double length;

	//read input file
	std::string junk;

	std::ifstream infile("input//INPUT_DIRECT_METHOD_FLOW.txt");
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
		return;
	}

	infile >> junk >> q >> junk >> k >> junk >> radius >> junk >> numseg;

	//print inputs to console
	std::cout << "Heat flow input: " << q << std::endl;
	std::cout << "Thermal conductivity: " << k << std::endl;
	std::cout << "Cylinder radius: " << radius << std::endl;
	std::cout << "Number of segments: " << numseg << std::endl;

	deltaT.resize(numseg, std::vector<double>(numseg));
	deltaU.resize(numseg, std::vector<double>(numseg));
	coordA.resize(2, std::vector<double>(numseg));
	coordB.resize(2, std::vector<double>(numseg));
	t0.resize(numseg);
	F.resize(numseg);
	coordPi.resize(2, std::vector<double>(numseg));

	//FILL THE COORDINATE VECTORS WITH THE SEGMENT START AND END POINTS
	//deltaTheta is the change in angle between segment normals
	deltaTheta = 2.0 * PI / numseg;
	//centerpoint of each segment is distance radius from the circle center
	//end points of each segment are distance radiusOuter from circle center since segments are straight. find radiusOuter with pythogoras
	radiusOuter = radius / cos(deltaTheta / 2.0);
	//let theta be a temp variable that gives angle to the starting point of each segment
	//loop over each segment, increment theta by deltaTheta and use to calculate the x and y coordinate of each segment start and end point
	//start with horizontal segment at top of circle. theta is measured from x-axis 
	double theta = (PI / 2.) - (deltaTheta / 2.0);
	//start point of segment 1
	coordA[0][0] = radiusOuter * cos(theta);
	coordA[1][0] = radiusOuter * sin(theta);
	for (size_t i = 0; i < numseg-1; i++) {
		theta += deltaTheta;
		coordB[0][i] = radiusOuter * cos(theta);
		coordB[1][i] = radiusOuter * sin(theta);
		coordA[0][i + 1] = coordB[0][i];
		coordA[1][i + 1] = coordB[1][i];
	}
	coordB[0][numseg - 1] = coordA[0][0];
	coordB[1][numseg - 1] = coordA[1][0];
	//calculate the coordinates of the Pi points at the center of each segment
	for (size_t i = 0; i < numseg; i++) {
		coordPi[0][i] = (coordA[0][i] + coordB[0][i]) / 2.0;
		coordPi[1][i] = (coordA[1][i] + coordB[1][i]) / 2.0;
	}

	//compute the applied tractions t0 due q at the center of element
	//start at element 1 (angle pi/2 from x-axis)
	theta = PI / 2.0;
	for (size_t i = 0; i < numseg; i++) {
		t0[i] = q * sin(theta);
		theta += deltaTheta;
	}

	//ASSEMBLE THE DELTA_T AND DELTA_U MATRICES
	//two loops. outer loop loops over each "source" segment Q
	//inner loop loops over each segment again, treating as the "potential" point Pi. loop over each segment and find the contribution of the source segment Q to each
	//better to do this way because can define the local coordinates fewer times
	
	double c = 0.5 / PI;
	double c1 = 0.5 / (PI * k);
	
	for (size_t i = 0; i < numseg; i++) {
		double dx = coordA[0][i] - coordB[0][i];
		double dy = coordA[1][i] - coordB[1][i];
		length = sqrt(dx * dx + dy * dy);
		//tangent vector
		Vt = { dx / length, dy / length };
		//normal vector
		Vn = { Vt[1], -1.0 * Vt[0] };
		for (size_t j = 0; j < numseg; j++) {
			//compute distance from Pi to start and end of segment rA and rB
			vrA[0] = coordA[0][i] - coordPi[0][j];
			vrA[1]  = coordA[1][i] - coordPi[1][j];
			vrB[0] = coordB[0][i] - coordPi[0][j];
			vrB[1] = coordB[1][i] - coordPi[1][j];
			double rA = sqrt(vrA[0] * vrA[0] + vrA[1] * vrA[1]);
			double rB = sqrt(vrB[0] * vrB[0] + vrB[1] * vrB[1]);
			//calculate the cosines and sines of the angles first, and then backcalculate the actual angle thetaA and thetaB
			//angles are in terms of local coordinates defined by the segment i. use the normal and tangent vectors to define the local x and y
			double cosThetaA = DotProduct(Vn, NormalizeVector(vrA, rA));
			double cosThetaB = DotProduct(Vn, NormalizeVector(vrB, rB));
			double sinThetaA = DotProduct(Vt, NormalizeVector(vrA, rA));
			double sinThetaB = DotProduct(Vt, NormalizeVector(vrB, rB));
			//compute angles
			double thetaA = acos(cosThetaA) * (sinThetaA >= 0 ? 1.0 : -1.0);
			double thetaB = acos(cosThetaB) * (sinThetaB >= 0 ? 1.0 : -1.0);

			//if diagonal entry treat differently than off-diagonal
			if (i == j) {
				deltaT[i][j] = 0.5;
				deltaU[i][j] = length * c1 * (log(length / 2.0) - 1.0);
			}
			else {
				deltaT[i][j] = c * (thetaB - thetaA);
				deltaU[i][j] = c1 * ((rB * sinThetaB * (log(rB) - 1) + thetaB * rB * cosThetaB) - (rA * sinThetaA * (log(rA) - 1) + thetaA * rA * cosThetaA));
			}
		}
	}

	//solve the system of equations by inverting the deltaT matrix and multiplying by F vector
	F = deltaU * t0;
	std::vector<std::vector<double>> invertedT = invertMatrix(deltaT);
	u = invertedT * F;

	//print boundary results to console
	std::cout << "Temperature on boundary segments: " << std::endl;
	for (size_t i = 0; i < numseg; i++) {
		std::cout << "Segment " << i << ": " << u[i]-q/k*coordPi[1][i] << std::endl;
	}
}
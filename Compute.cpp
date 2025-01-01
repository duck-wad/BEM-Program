#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

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

	double q, k, radius, radiusOuter, deltaTheta;
	int numseg;
	//kernel matrices
	std::vector<std::vector<double>> deltaT;
	std::vector<std::vector<double>> deltaU;
	//vector with segment center potential
	Eigen::VectorXd u;
	//vector for deltaU times t0 (applied flow)
	std::vector<double> F;
	//applied flow i.e boundary condition
	std::vector<double> t0;
	//start coordA and end coordB coordinates of each segment
	//first row contains the x coord, second row contains the y coord
	//each column correspond with a point
	std::vector<std::vector<double>> coordA, coordB;
	//coordinate of the center point Pi
	//used for solution and postprocessing
	std::vector<std::vector<double>> coordPi;
	//temp tangent vector for each segment
	std::vector<std::vector<double>> Vt;
	//temp normal unit vector for each segment
	std::vector<std::vector<double>> Vn;
	//temp vectors pointing from Pi to A and B of segment
	//also used for postprocessing
	std::vector<double> vrA(2);
	std::vector<double> vrB(2);
	//segment length
	double length;
	//number of domain points to compute temp and flow
	int inPoints = 0;
	//coordinates of interior points, first row x second row y
	std::vector<std::vector<double>> coordPa;
	//vector containing domain point temperatures
	std::vector<double> uP;
	//vector containing domain point flows, first row x second row y
	std::vector<std::vector<double>> qP;

	//read input file
	std::string junk;

	std::ifstream infile("input//INPUT_DIRECT_METHOD_FLOW.txt");
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
		return;
	}

	infile >> junk >> q >> junk >> k >> junk >> radius >> junk >> numseg;
	infile >> junk;
	if (junk == "inpoint") {
		infile >> junk >> inPoints;
		coordPa.resize(2, std::vector<double>(inPoints));
		for (size_t i = 0; i < inPoints; i++) {
			infile >> junk >> coordPa[0][i] >> coordPa[1][i];
		}
	}

	infile.close();

	deltaT.resize(numseg, std::vector<double>(numseg));
	deltaU.resize(numseg, std::vector<double>(numseg));
	coordA.resize(2, std::vector<double>(numseg));
	coordB.resize(2, std::vector<double>(numseg));
	t0.resize(numseg);
	F.resize(numseg);
	coordPi.resize(2, std::vector<double>(numseg));
	//different form than coord vectors because need to use DotProduct function
	Vt.resize(numseg, std::vector<double>(2));
	Vn.resize(numseg, std::vector<double>(2));

	/* FILL THE COORDINATE VECTORS WITH THE SEGMENT START AND END POINTS */

	//deltaTheta is the change in angle between segment normals
	deltaTheta = 2.0 * PI / numseg;
	//centerpoint of each segment is distance radius from the circle center
	//end points of each segment are distance radiusOuter from circle center since segments are straight. find radiusOuter with pythogoras
	radiusOuter = radius / cos(deltaTheta / 2.0);
	//let theta be a temp variable that gives angle to the starting point of each segment
	//loop over each segment, increment theta by deltaTheta and use to calculate the x and y coordinate of each segment start and end point
	//start with horizontal segment at top of circle. theta is measured from x-axis 
	double theta = (PI / 2.0) - (deltaTheta / 2.0);
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

	/* ASSEMBLE THE DELTA_T AND DELTA_U MATRICES */

	//two loops. outer loop loops over each "source" segment Q
	//inner loop loops over each segment again, treating as the "potential" point Pi. loop over each segment and find the contribution of the source segment Q to each
	double c = 0.5 / PI;
	double c1 = 0.5 / (PI * k);
	
	for (size_t i = 0; i < numseg; i++) {
		double dx = coordA[0][i] - coordB[0][i];
		double dy = coordA[1][i] - coordB[1][i];
		length = sqrt(dx * dx + dy * dy);
		//tangent vector
		Vt[i] = {(dx / length), (dy / length)};
		//normal vector
		Vn[i] = {Vt[i][1], (-1.0 * Vt[i][0])};
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
			double cosThetaA = DotProduct(Vn[i], vrA) / rA;
			double cosThetaB = DotProduct(Vn[i], vrB) / rB;
			double sinThetaA = DotProduct(Vt[i], vrA) / rA;
			double sinThetaB = DotProduct(Vt[i], vrB) / rB;
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

	F = deltaU * t0;
	
	/* SOLVE SYSTEM USING EIGEN LIBRARY */

	// Convert to Eigen matrices and vectors
	Eigen::MatrixXd deltaT_eigen(numseg, numseg);
	Eigen::MatrixXd deltaU_eigen(numseg, numseg);
	Eigen::VectorXd F_eigen(numseg);

	for (size_t i = 0; i < numseg; ++i) {
		F_eigen(i) = F[i];
		for (size_t j = 0; j < numseg; ++j) {
			deltaT_eigen(i, j) = deltaT[i][j];
			deltaU_eigen(i, j) = deltaU[i][j];
		}
	}

	// Regularization, since the system is sensitive to input precision
	double lambda = 1e-6;  // Small regularization parameter
	Eigen::MatrixXd deltaT_regularized = deltaT_eigen + lambda * Eigen::MatrixXd::Identity(deltaT_eigen.rows(), deltaT_eigen.cols());

	// Solve with regularization
	u = deltaT_regularized.colPivHouseholderQr().solve(F_eigen);

	/* COMPUTE RESULTS IN DOMAIN (OUTSIDE OF CYLINDER) */

	uP.resize(inPoints, 0.0);
	qP.resize(2, std::vector<double>(inPoints, 0.0));

	if (inPoints > 0) {
		for (size_t i = 0; i < inPoints; i++) {
			//for each domain point, need to loop over each segment of boundary to get contribution from each segment to flow and potential 
			for (size_t j = 0; j < numseg; j++) {
				//need angles between the boundary point and domain point
				//calculated same way as previous but the sign on the normal vector is times -1 since normal vec points inwards, and domain is outside of cylinder
				vrA[0] = coordA[0][j] - coordPa[0][i];
				vrA[1] = coordA[1][j] - coordPa[1][i];
				vrB[0] = coordB[0][j] - coordPa[0][i];
				vrB[1] = coordB[1][j] - coordPa[1][i];
				double rA = sqrt(vrA[0] * vrA[0] + vrA[1] * vrA[1]);
				double rB = sqrt(vrB[0] * vrB[0] + vrB[1] * vrB[1]);

				double cosThetaA = -1 * DotProduct(Vn[j], vrA) / rA;
				double cosThetaB = -1 * DotProduct(Vn[j], vrB) / rB;
				double sinThetaA = -1 * DotProduct(Vt[j], vrA) / rA;
				double sinThetaB = -1 * DotProduct(Vt[j], vrB) / rB;
				double h = rA * cosThetaA;
				//compute angles
				double thetaA = acos(cosThetaA) * (sinThetaA >= 0 ? 1.0 : -1.0);
				double thetaB = acos(cosThetaB) * (sinThetaB >= 0 ? 1.0 : -1.0);
				//adjust thetaA to have thetaB-thetaA lower than 180deg
				if ((thetaB - thetaA) > PI) {
					thetaA += 2.0 * PI;
				}

				//compute the domain point temperature by adding contributions from each segment
				//use deltaT, deltaU, and segment temp and tractions to find potential. definition of dT and dU same as before
				double dT = c * (thetaB - thetaA);
				double dU = c1 * ((rB * sinThetaB * (log(rB) - 1) + thetaB * rB * cosThetaB) - (rA * sinThetaA * (log(rA) - 1) + thetaA * rA * cosThetaA));
				uP[i] += dU * t0[j] - dT * u[j];

				//compute the domain point flows in x and y direction
				double dSx = c1 * (thetaB - thetaA);
				double dSy = 0.0;
				double frac = cosThetaB / cosThetaA;
				if (frac > 0.0) {
					dSy = -1 * c1 * log(frac);
				}
				double dRx = -1 * c / h * (cosThetaB * sinThetaB - cosThetaA * sinThetaA);
				double dRy = c / h * (cosThetaB * cosThetaB - cosThetaA * cosThetaA);
				//compute the flows in terms of local coordinates first
				double qPx = -1 * k * (dSx * t0[j] - dRx * u[j]);
				double qPy = -1 * k * (dSy * t0[j] - dRy * u[j]);
				//convert to global coordinates and add the contributions
				qP[0][i] += qPx * Vn[j][0] - qPy * Vn[j][1];
				qP[1][i] += qPx * Vn[j][1] + qPy * Vn[j][0];
			}
			//make some sort of adjustment idek why. flow in y need to add the applied flow to the domain from the input file
			uP[i] -= q / k * coordPa[1][i];
			qP[1][i] += q;
		}
	}

	/* PRINT RESULTS */

	std::cout << "Heat flow past cylinder using the direct BE method." << std::endl;

	//print inputs to console
	std::cout << "Heat flow input: " << q << std::endl;
	std::cout << "Thermal conductivity: " << k << std::endl;
	std::cout << "Cylinder radius: " << radius << std::endl;
	std::cout << "Number of segments: " << numseg << std::endl;

	//print boundary results to console
	std::cout << "Temperature at segment centers: " << std::endl;
	for (size_t i = 0; i < numseg; i++) {
		std::cout << "Segment " << i + 1 << ": " << u[i] - (q / k * coordPi[1][i]) << std::endl;
	}

	//print domain results to console
	std::cout << "Temperature and flow at domain points: " << std::endl;
	for (size_t i = 0; i < inPoints; i++) {
		std::cout << "Coordinates: x=" << coordPa[0][i] << " y=" << coordPa[1][i] << "   ";
		std::cout << "T=" << uP[i] << "   ";
		std::cout << "q-x=" << qP[0][i] << "   q-y=" << qP[1][i];
		std::cout << std::endl;
	}
}
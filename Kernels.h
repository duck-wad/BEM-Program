#pragma once
#include <vector>

//computes fundamental solution to Laplace problems for potential and flow
double Potential(double r, double k, int Cdim);
double Flow(double r, const std::vector<double>& dxr, const std::vector<double>& normal, int Cdim);

//computes the fundamental solution for displacements and tractions (Kelvin solution) for the 2D plane strain and 3D case
void KelvinDisplacement(std::vector<std::vector<double>>& UK, const std::vector<double>& dxr, double r, double E, double nu, int Cdim);
void KelvinTraction(std::vector<std::vector<double>>& TK, const std::vector<double>& dxr, const std::vector<double>& normal, double r, double nu, int Cdim);
#include <iostream>

#include "Utils.h"

std::vector<std::vector<double>>& operator*=(std::vector<std::vector<double>>& matrix, double constant) {

	for (auto& row : matrix) {
		for (auto& value : row) {
			value *= constant;
		}
	}
	return matrix;
}

void PrintMatrix(const std::vector<std::vector<double>>& matrix){
	for (const auto& row : matrix) {
		for (const auto& value : row) {
			std::cout << value << "\t";
		}
		std::cout << std::endl;
	}
}
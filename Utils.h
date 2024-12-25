#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

const double PI = 3.14159265359;

const double LOW_TOL = 1e-8;
const double HIGH_TOL = 1e-12;

int Kronecker(int i, int j);

template<typename T>
void PrintMatrix(const std::vector<std::vector<T>>& matrix) {
	for (const auto& row : matrix) {
		for (const auto& value : row) {
			std::cout << value << "\t";
		}
		std::cout << std::endl;
	}
}

template<typename T>
void writeVectorToCSV(const std::vector<T>& vector, const std::string& filename) {

	std::ofstream file(filename);

	if (!file.is_open()) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	for (size_t i = 0; i < vector.size(); i++) {
		file << vector[i] << "\n";
	}

	file.close();

	std::cout << "Vector written to " << filename << " succesfully." << std::endl;
}

template<typename T>
void writeMatrixToCSV(const std::vector<std::vector<T>>& matrix, const std::string& filename) {

	// Open the file stream

	std::ofstream file(filename);

	if (!file.is_open()) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	size_t rows = matrix[0].size();
	size_t cols = matrix.size();

	// Write the matrix to the file in row-major order
	for (size_t i = 0; i < cols; ++i) {
		for (size_t j = 0; j < rows; ++j) {
			file << matrix[i][j];
			if (j < rows - 1) { // Add a comma unless it's the last column
				file << ",";
			}
		}
		file << "\n"; // Newline after each row
	}

	file.close();

	std::cout << "Matrix written to " << filename << " successfully." << std::endl;
}

// Function to read a CSV file into a std::vector<std::vector<T>>
template <typename T>
std::vector<std::vector<T>> readCSVtoMatrix(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	std::vector<std::vector<T>> data;
	std::string line;
	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<T> row;
		while (std::getline(lineStream, cell, ',')) {
			row.push_back(static_cast<T>(std::stod(cell)));
		}
		data.push_back(row);
	}
	return data;
}

// Function to read a CSV file into a std::vector<T>
template <typename T>
std::vector<T> readCSVtoVector(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	std::vector<T> data;
	std::string line;
	while (std::getline(file, line)) {
		data.push_back(static_cast<T>(std::stod(line)));
	}
	return data;
}

/* ALGEBRAIC MATRIX / VECTOR OPERATOR OVERLOADS */

template<typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}

	size_t mat1row = mat1.size();
	size_t mat1col = mat1[0].size();
	size_t mat2row = mat2.size();
	size_t mat2col = mat2[0].size();

	if (mat1col != mat2row || mat1row != mat2col) {
		throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
	}

	std::vector<std::vector<T>> output(mat2col, std::vector<T>(mat1row, T()));

	for (size_t i = 0; i < mat2col; i++) {
		for (size_t j = 0; j < mat1row; j++) {
			for (size_t k = 0; k < mat1col; k++) {
				output[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
	return output;
}

template<typename T>
std::vector<T> operator* (const std::vector<std::vector<T>>& mat, const std::vector<T>& vec) {
	if (mat.empty() || vec.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat.size() != vec.size() || mat[0].size() != vec.size()) {
		throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
	}
	std::vector<T> output(vec.size(), T());

	for (size_t i = 0; i < vec.size(); i++) {
		for (size_t j = 0; j < vec.size(); j++) {
			output[i] += mat[i][j] * vec[j];
		}
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& mat, const T c) {
	if (mat.empty()) {
		throw std::invalid_argument("Matrix cannot be empty");
	}

	std::vector<std::vector<T>> output(mat.size(), std::vector<T>(mat[0].size(), 0.0));

	for (size_t i = 0; i < mat.size(); i++) {
		for (size_t j = 0; j < mat[i].size(); j++) {
			output[i][j] = mat[i][j] * c;
		}
	}
	return output;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T c) {
	if (vec.empty()) {
		throw std::invalid_argument("Vector cannot be empty");
	}
	std::vector<T> output(vec.size(), 0.0);
	for (size_t i = 0; i < vec.size(); i++) {
		output[i] += vec[i] * c;
	}
	return output;
}

template<typename T>
std::vector<T> operator/(const std::vector<T>& vec, const T c) {
	if (vec.empty()) {
		throw std::invalid_argument("Vector cannot be empty");
	}
	std::vector<T> output(vec.size());
	for (size_t i = 0; i < vec.size(); i++) {
		output[i] += vec[i] / c;
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>>& operator*= (std::vector<std::vector<T>>& matrix, T scalar) {
	for (auto& row : matrix) {
		for (auto& value : row) {
			value *= scalar;
		}
	}
	return matrix;
}

template<typename T>
std::vector<std::vector<T>>& operator/= (std::vector<std::vector<T>>& matrix, T scalar) {
	for (auto& row : matrix) {
		for (auto& value : row) {
			value /= scalar;
		}
	}
	return matrix;
}

template<typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	std::vector<std::vector<T>> output(mat1.size(), std::vector<T>(mat1[0].size(), T()));

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			output[i][j] = mat1[i][j] + mat2[i][j];
		}
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	std::vector<std::vector<T>> output(mat1.size(), std::vector<T>(mat1[0].size(), T()));

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			output[i][j] = mat1[i][j] - mat2[i][j];
		}
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			mat1[i][j] += mat2[i][j];
		}
	}
	return mat1;
}

template<typename T>
std::vector<std::vector<T>>& operator-=(std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			mat1[i][j] -= mat2[i][j];
		}
	}
	return mat1;
}

template<typename T>
std::vector<T>& operator+=(std::vector<T>& v1, const std::vector<T>& v2) {

	if (v1.empty() || v2.empty()) {
		throw std::invalid_argument("Vectors must not be empty");
	}
	if (v1.size() != v2.size()) {
		throw std::invalid_argument("Vectors must be the same size for addition");
	}

	for (size_t i = 0; i < v1.size(); i++) {
		v1[i] += v2[i];
	}
	return v1;
}

template<typename T>
std::vector<T>& operator-=(std::vector<T>& v1, const std::vector<T>& v2) {

	if (v1.empty() || v2.empty()) {
		throw std::invalid_argument("Vectors must not be empty");
	}
	if (v1.size() != v2.size()) {
		throw std::invalid_argument("Vectors must be the same size for addition");
	}

	for (size_t i = 0; i < v1.size(); i++) {
		v1[i] -= v2[i];
	}
	return v1;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {

	if (v1.empty() || v2.empty()) {
		throw std::invalid_argument("Vectors must not be empty");
	}
	if (v1.size() != v2.size()) {
		throw std::invalid_argument("Vectors must be the same size for addition");
	}

	std::vector<T> output(v1.size());

	for (size_t i = 0; i < v1.size(); i++) {
		output[i] = v1[i] + v2[i];
	}
	return output;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2) {

	if (v1.empty() || v2.empty()) {
		throw std::invalid_argument("Vectors must not be empty");
	}
	if (v1.size() != v2.size()) {
		throw std::invalid_argument("Vectors must be the same size for addition");
	}

	std::vector<T> output(v1.size());

	for (size_t i = 0; i < v1.size(); i++) {
		output[i] = v1[i] - v2[i];
	}
	return output;
}

/* END OF ALGEBRAIC MATRIX OPERATOR OVERLOADS */
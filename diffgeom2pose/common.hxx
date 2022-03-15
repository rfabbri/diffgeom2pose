/* Includes common functions to all modules */
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

namespace common {
	// TODO: check if `colon()` is correct from MATLAB documentation
	template<typename T>
	T *colon(T start, T step, T end)
	{
		int num_elements = floor((end - start) / step);
		static T output[num_elements];

		for (int i = 0; i < num_elements; i++) {
			output[i] = start * i + step;
		}

		return output;
	}

	// TODO: optimize! (make inline or use #define)
	template<typename T>
	T norm(T *input_vector, int length)
	{
		T output;
		for (int i = 0; i < length; i++) {
			output += abs(input[i] * input[i]);
		}
		return output;
	}

	template<typename T>
	T det3x3(T* input_matrix, int rows, int cols)
	{
		// TODO: Check if this algorithm is correct
		// Matrix must be 3x3 square
		if (rows != cols != 3) {
			return nullptr;
		}
		T x = ((input_matrix[1][1] * input_matrix[2][2]) - (input_matrix[2][1] * input_matrix[1][2]));
		T y = ((input_matrix[1][0] * input_matrix[2][2]) - (input_matrix[2][0] * input_matrix[1][2]));
		T z = ((input_matrix[1][0] * input_matrix[2][1]) - (input_matrix[2][0] * input_matrix[1][1]));

		return input_matrix[0][0] * x - input_matrix[0][1] * y - input_matrix[0][2] * z;

	}

	template<typename T>
	T intpow(T base, int exponent)
	{
		T output = 1;
		while (exponent-- > 0) {
			output *= base;
		}
		return output;
	}

	template<typename T>
	T **expm(T **M) 
	{
		// TODO: implement `expm()` (matrix exponential)
	}

	// TODO: implement exception case for single column or single row
	// return as vector
	template<typename T>
	T **randn(int rows, int cols) 
	{
		// TODO: check if normal distribution implementation is correct
		T **output = (T**)malloc(rows * sizeof(T*));
		for (int i = 0; i < rows; i++) {
			output[i] = (T*)malloc(cols * sizeof(T));
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<T> d(0,1);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				output[i][j] = d(gen);
			}
		}

		return output;
	}

	template<typename T>
	T *vec_subtract(T vec1[], int vec1_size, T vec2[], int vec2_size)
	{
		// Vectors must have same size
		if (vec1_size != vec2_size) {
			return nullptr;
		}

		static T output[vec1_size];

		for (int i = 0; i < vec1_size; i++) {
			output[i] = vec1[i] - vec2[i];
		}

		return output;
	}
}


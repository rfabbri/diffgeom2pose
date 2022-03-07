/* Includes common functions to all modules */
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

namespace common {
	// TODO: check if `colon()` is correct from MATLAB documentation
	template<typename T>
	T* colon(T start, T step, T end)
	{
		int num_elements = floor((end - start) / step);
		T* output = (T*)malloc(sizeof(T * num_elements));

		for (int i = 0; i < num_elements; i++) {
			output[i] = start * i + step;
		}

		return output;
	}

	// TODO: optimize! (make inline or use #define)
	template<typename T>
	T norm(T* input_vector, int length)
	{
		T output;
		for (int i = 0; i < length; i++) {
			output += abs(input[i] * input[i]);
		}
		return output;
	}

	template<typename T>
	T det(T* input_matrix) 
	{
		// TODO: implement `det()`	
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
	T** expm(T** M) 
	{
		// TODO: implement `expm()` (matrix exponential)
	}

	// TODO: implement exception case for single column or single row
	// return as vector
	template<typename T>
	T** randn(int rows, int cols) 
	{
		// TODO: check if normal distribution implementation is correct
		T** output = (T**)malloc(rows * sizeof(T*));
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
}


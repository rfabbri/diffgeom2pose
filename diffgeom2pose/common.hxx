/* Includes common functions to all modules */
#include <iostream>
#include <cmath>

namespace common {
	// TODO: check if `colon()` is correct from MATLAB documentation
	template<typename T>
	void colon(T start, T step, T end, T output[])
	{
		int num_elements = floor((end - start) / step);

		for (int i = 0; i < num_elements; i++) {
			output[i] = start * i + step;
		}

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
	T det3x3(T input_matrix[3][3])
	{
		// TODO: Check if this algorithm is correct
		// Matrix must be 3x3 square
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

	template<typename T>
	void vec_sub(T vec1[3], T vec2[3], T output[3])
	{
		// For vectors of 3 elements
		for (int i = 0; i < 3; i++) {
			output[i] = vec1[i] - vec2[i];
		}
	}

	template<typename T>
	void vec_el_wise_pow(T vec[2001], T exp, T output[2001])
	{
		for (int i = 0; i < 2001; i++) {
			output[i] = pow(vec[i], exp);
		}
	}

	template<typename T>
	void vec_el_wise_mult(T vec1[3], T vec2[3], T output[3])
	{
		for (int i = 0; i < 3; i++) {
			output[i] = vec1[i] * vec2[i];
		}
	}

	template<typename T>
	void vec_mult_by_scalar(T vec[2001], T scalar, T output[2001])
	{
		for (int i = 0; i < 2001; i++) {
			output[i] = vec[i] * scalar;
		}
	}

	template<typename T>
	T vec_sum(T vec[3])
	{
		T output = 0;
		for (int i = 0; i < 3; i++) {
			output += vec[i];
		}
		return output;
	}

}


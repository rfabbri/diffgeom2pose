/* Includes common functions to all modules */
#include <iostream>
#include <vector>
#include <cmath>

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
T int_pow(T base, int exponent)
{
	T output = 1;
	for (int i = 0; i < exponent; i++) {
		output *= base;
	}
	return output;
}

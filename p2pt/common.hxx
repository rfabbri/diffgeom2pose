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
	T norm(T *input_vec, int length)
	{
		T output = 0;
		for (int i = 0; i < length; i++) {
			output += abs(input_vec[i] * input_vec[i]);
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

	// TODO: Make Matrix into a class; simplify several operations
	template<typename T>
	void invm3x3(T input_m[3][3], T output_m[3][3])
	{
		// 3x3 MATRIX INVERSION ALGORITHM
		//             -1               T
		//  -1  [a b c]      1   [A B C]      1   [A D G]
		// A  = [d e f] = ------ [D E F] = ------ [B E H]
		//      [g h i]   det(A) [G H I]   det(A) [C F I]
		//
		// A =  (ei - fh), D = -(bi - ch), G =  (bf - ce),
		// B = -(di - fg), E =  (ai - cg), H = -(af - cd),
		// C =  (dh - eg), F = -(ah - bg), I =  (ae - bd).

		static T
		a_ = input_m[0][0], b_ = input_m[0][1], c_ = input_m[0][2],
		d_ = input_m[1][0], e_ = input_m[1][1], f_ = input_m[1][2],
		g_ = input_m[2][0], h_ = input_m[2][1], i_ = input_m[2][2];

		static T
		A_ =  (e_*i_ - f_*h_), B_ = -(d_*i_ - f_*g_), C_ =  (d_*h_ - e_*g_),
		D_ = -(b_*i_ - c_*h_), E_ =  (a_*i_ - c_*g_), F_ = -(a_*h_ - b_*g_),
		G_ =  (b_*f_ - c_*e_), H_ = -(a_*f_ - c_*d_), I_ =  (a_*e_ - b_*d_);

		static T invdet_A = 1 / (a_*A_ + b_*B_ + c_*C_);

		output_m[0][0] = invdet_A * A_; output_m[0][1] = invdet_A * D_; output_m[0][2] = invdet_A * G_;
		output_m[1][0] = invdet_A * B_; output_m[1][1] = invdet_A * E_; output_m[1][2] = invdet_A * H_;
		output_m[2][0] = invdet_A * C_; output_m[2][1] = invdet_A * F_; output_m[2][2] = invdet_A * I_;
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
	T vec1vec2_dot(T vec1[3], T vec2[3])
	{
		T output = 0;
		for (int i = 0; i < 3; i++) {
			output += vec1[i] * vec2[i];
		}
		return output;
	}

	template<typename T>
	void vec1vec2_sum(T vec1[t_vector_len], T vec2[t_vector_len], T output[t_vector_len])
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] + vec2[i];
		}
	}

	template<typename T>
	void vec1vec2_sub(const T vec1[3], const T vec2[3], T output[3])
	{
		// For vectors of 3 elements
		for (int i = 0; i < 3; i++) {
			output[i] = vec1[i] - vec2[i];
		}
	}

	template<typename T>
	void vec_add_scalar(T vec1[t_vector_len], int scalar, T output[t_vector_len])
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] + scalar;
		}
	}

	template<typename T>
	void vec_el_wise_pow(const T vec[t_vector_len], int exp, T output[t_vector_len])
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = pow(vec[i], exp);
		}
	}

	template<typename T>
	void vec1vec2_el_wise_right_div(T vec1[t_vector_len], T vec2[t_vector_len], T output[t_vector_len])
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] / vec2[i];
		}
	}


	// Returns pointer for simplicity of calculations in `rf_sample_pose_poly`
	// TODO: Check if there is a better way to do this (variadic function perhaps?)
	#pragma region vec_el_wise_mult

	// Special case for 2 vectors of 3 elements each
	// TODO: Check if there is a batter way to do this
	template<typename T>
	T* vec_3el_wise_mult2(
		const T vec1[t_vector_len], const T vec2[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < 3; i++) {
			output[i] = vec1[i] * vec2[i];
		}
		return output;
	}

	// For 3 vectors
	template<typename T>
	T *vec_el_wise_mult3(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i];
		}
		return output;
	}

	// For 4 vectors
	template<typename T>
	T *vec_el_wise_mult4(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len], const T vec4[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i] * vec4[i];
		}
		return output;
	}

	// For 5 vectors
	template<typename T>
	T *vec_el_wise_mult5(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len], const T vec4[t_vector_len],
		const T vec5[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i] * vec4[i]
				* vec5[i];
		}
		return output;
	}

	// For 6 vectors
	template<typename T>
	T *vec_el_wise_mult6(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len], const T vec4[t_vector_len],
		const T vec5[t_vector_len], const T vec6[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i] * vec4[i]
				* vec5[i] * vec6[i];
		}
		return output;
	}

	// For 7 vectors
	template<typename T>
	T *vec_el_wise_mult7(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len], const T vec4[t_vector_len],
		const T vec5[t_vector_len], const T vec6[t_vector_len], const T vec7[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i] * vec4[i]
				* vec5[i] * vec6[i] * vec7[i];
		}
		return output;
	}

	// For 8 vectors
	template<typename T>
	T *vec_el_wise_mult8(
		const T vec1[t_vector_len], const T vec2[t_vector_len], const T vec3[t_vector_len], const T vec4[t_vector_len],
		const T vec5[t_vector_len], const T vec6[t_vector_len], const T vec7[t_vector_len], const T vec8[t_vector_len],
		T output[t_vector_len]
	)
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec1[i] * vec2[i] * vec3[i] * vec4[i]
				* vec5[i] * vec6[i] * vec7[i] * vec8[i];
		}
		return output;
	}
	#pragma endregion

	template<typename T>
	void vec_mult_by_scalar(int scalar, T vec[t_vector_len], T output[t_vector_len])
	{
		for (int i = 0; i < t_vector_len; i++) {
			output[i] = vec[i] * scalar;
		}
	}

	template<typename T>
	T vec_el_sum(T vec[3])
	{
		T output = 0;
		for (int i = 0; i < 3; i++) {
			output += vec[i];
		}
		return output;
	}

}


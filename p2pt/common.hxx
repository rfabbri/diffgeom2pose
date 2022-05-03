/* Includes common functions to all modules */
#ifndef common_hxx_
#define common_hxx_

#include <cmath>

namespace common {

/* unused
template<typename T>
inline void
colon(T start, T step, T end, T output[])
{
	int num_elements = floor((end - start) / step) + 1;

	for (int i = 0; i < num_elements; i++)
		output[i] = start + i * step;

}

template<typename T>
inline T
norm(const T *input_vec, int length)
{
	T output = 0;
	for (int i = 0; i < length; i++) {
		output += std::abs(input_vec[i] * input_vec[i]);
	}
	return std::sqrt(output);
}
*/

template<typename T>
inline T
det3x3(const T (&m)[3][3])
{
	return (m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1])
		 - (m[2][0]*m[1][1]*m[0][2] + m[2][1]*m[1][2]*m[0][0] + m[2][2]*m[1][0]*m[0][1]);
}

template<typename T>
inline void
invm3x3(const T (&input_m)[3][3], T (&output_m)[3][3])
{
	// 3x3 MATRIX INVERSION ALGORITHM
	//             -1               T
	//  -1  [a b c]      1   [A B C]      1   [A D G]
	// M  = [d e f] = ------ [D E F] = ------ [B E H]
	//      [g h i]   det(M) [G H I]   det(M) [C F I]
	//
	// A =  (ei - fh), D = -(bi - ch), G =  (bf - ce),
	// B = -(di - fg), E =  (ai - cg), H = -(af - cd),
	// C =  (dh - eg), F = -(ah - bg), I =  (ae - bd).

	const T 
	a = input_m[0][0], b = input_m[0][1], c = input_m[0][2],
	d = input_m[1][0], e = input_m[1][1], f = input_m[1][2],
	g = input_m[2][0], h = input_m[2][1], i = input_m[2][2];

	const T 
	A =  (e*i - f*h), B = -(d*i - f*g), C =  (d*h - e*g),
	D = -(b*i - c*h), E =  (a*i - c*g), F = -(a*h - b*g),
	G =  (b*f - c*e), H = -(a*f - c*d), I =  (a*e - b*d);

	const T invdet_M = 1. / (a*A + b*B + c*C);

	output_m[0][0] = invdet_M * A; output_m[0][1] = invdet_M * D; output_m[0][2] = invdet_M * G;
	output_m[1][0] = invdet_M * B; output_m[1][1] = invdet_M * E; output_m[1][2] = invdet_M * H;
	output_m[2][0] = invdet_M * C; output_m[2][1] = invdet_M * F; output_m[2][2] = invdet_M * I;
}

template<typename T>
inline void
multm3x3(const T (&m1)[3][3], const T (&m2)[3][3], T output_m[][3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			output_m[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				output_m[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

template<typename T>
inline void
multm_3x3_3x1(const T m1[][3], const T (&m2)[3], T (&output_m)[3])
{
	for (int i = 0; i < 3; i++) {
		output_m[i] = 0;
		for (int j = 0; j < 3; j++) {
			output_m[i] += m1[i][j] * m2[j];
		}
	}
}

template<typename T>
inline T
vec1vec2_3el_dot(const T (&vec1)[3], const T (&vec2)[3])
{
	T output = 0;
	for (int i = 0; i < 3; i++)
		output += vec1[i] * vec2[i];
	return output;
}

/* not used
template<typename T>
inline void
vec1vec2_sum(T (&vec1)[T_VECTOR_LEN], T (&vec2)[T_VECTOR_LEN], T (&output)[T_VECTOR_LEN])
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] = vec1[i] + vec2[i];
}
*/

template<typename T>
inline void
vec1vec2_3el_sum(T (&vec1)[3], T (&vec2)[3], T (&output)[3])
{
		output[0] = vec1[0] + vec2[0];
		output[1] = vec1[1] + vec2[1];
		output[2] = vec1[2] + vec2[2];
}

template<typename T>
inline void
vec1vec2_3el_sub(const T (&vec1)[3], const T (&vec2)[3], T (&output)[3])
{
	// For vectors of 3 elements
		output[0] = vec1[0] - vec2[0];
		output[1] = vec1[1] - vec2[1];
		output[2] = vec1[2] - vec2[2];
}

/* not used
template<typename T>
inline void
vec_add_scalar(T (&vec1)[T_VECTOR_LEN], int scalar, T (&output)[T_VECTOR_LEN])
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] = vec1[i] + scalar;
}

template<typename T>
inline void
vec_el_wise_pow(const T (&__restrict__ vec)[T_VECTOR_LEN], int exp, T (&__restrict__ output)[T_VECTOR_LEN])
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] = intpow(vec[i], exp);
}

template<typename T>
inline void
vec1vec2_el_wise_right_div(T (&vec1)[T_VECTOR_LEN], T (&vec2)[T_VECTOR_LEN], T (&output)[T_VECTOR_LEN])
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] = vec1[i] / vec2[i];
}
*/

template<typename T>
inline void
vec_3el_div_by_scalar(const T scalar, T (&vec)[3], T (&output)[3])
{
	for (int i = 0; i < 3; i++)
		output[i] = vec[i] / scalar;
}

// Special case for 2 vectors of 3 elements each
// Returns pointer for simplicity of calculations in `sample_pose_poly`
template<typename T>
inline T
*vec_3el_wise_mult2(
	const T (&vec1)[3], const T (&vec2)[3],
	T (&output)[3]
)
{
	output[0] = vec1[0] * vec2[0];
	output[1] = vec1[1] * vec2[1];
	output[2] = vec1[2] * vec2[2];
	return output;
}

// For 3 vectors and add
/* not used 
template<typename T>
inline void
vec_el_wise_mult3_and_add(
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i];
}

// For 4 vectors and add
template<typename T>
inline void
vec_el_wise_mult4_and_add(
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i];
}

// For 5 vectors multiply by scalar and add
template<typename T>
inline void
vec_el_wise_mult5_scalar_and_add(
	const int scalar,
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i] * scalar;
}

// For 5 vectors and add
template<typename T>
inline void
vec_el_wise_mult5_and_add(
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i];
}

// For 6 vectors multiply by scalar and add
template<typename T>
inline void
vec_el_wise_mult6_scalar_and_add(
	const int scalar,
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN], const T (&__restrict__ vec6)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i] * vec6[i] * scalar;
}

// For 6 vectors and add
template<typename T>
inline void
vec_el_wise_mult6_and_add(
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN], const T (&__restrict__ vec6)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i] * vec6[i];
}

// For 7 vectors multiply by scalar and add
template<typename T>
inline void
vec_el_wise_mult7_scalar_and_add(
	const int scalar,
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN], const T (&__restrict__ vec6)[T_VECTOR_LEN], const T (&__restrict__ vec7)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i] * vec6[i] * vec7[i] * scalar;
}

// For 8 vectors multiply by scalar and add
template<typename T>
inline void
vec_el_wise_mult8_scalar_and_add(
	const int scalar,
	const T (&__restrict__ vec1)[T_VECTOR_LEN], const T (&__restrict__ vec2)[T_VECTOR_LEN], const T (&__restrict__ vec3)[T_VECTOR_LEN], const T (&__restrict__ vec4)[T_VECTOR_LEN],
	const T (&__restrict__ vec5)[T_VECTOR_LEN], const T (&__restrict__ vec6)[T_VECTOR_LEN], const T (&__restrict__ vec7)[T_VECTOR_LEN], const T (&__restrict__ vec8)[T_VECTOR_LEN],
	T (&__restrict__ output)[T_VECTOR_LEN]
)
{
	for (int i = 0; i < T_VECTOR_LEN; i++)
		output[i] += vec1[i] * vec2[i] * vec3[i] * vec4[i]
		           * vec5[i] * vec6[i] * vec7[i] * vec8[i] * scalar;
}
*/

template<typename T>
inline void
vec_3el_mult_by_scalar(const T scalar, const T (&vec)[3], T (&output)[3])
{
	for (int i = 0; i < 3; i++)
		output[i] = vec[i] * scalar;
}

template<typename T>
inline T
vec_3el_sum(T (&vec)[3])
{
	T output = 0;
	for (int i = 0; i < 3; i++)
		output += vec[i];
	return output;
}

}

#endif // common_hxx_


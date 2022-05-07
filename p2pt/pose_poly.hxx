#ifndef pose_poly_hxx_
#define pose_poly_hxx_

#include "pose_poly.h"

#include <complex>

namespace P2Pt {
  
template<typename T>
void
pose_poly<T>::
get_sigmas(const int ts_len, const T (&ts)[ROOT_IDS_LEN],
	T (*out)[2][TS_MAX_LEN][TS_MAX_LEN], int (*out_len)[2][TS_MAX_LEN])
{
	/*
	 `out`
	       [0] -> sigmas1[TS_MAX_LEN][TS_MAX_LEN]
	       [1] -> sigmas2[TS_MAX_LEN][TS_MAX_LEN]

	       `sigmasX` (can contain single values or array of values)
	             [0][0] -> float/double
	             [1][0] -> float/double
	             [2][0] -> float/double
	             [3][ ] -> [0] = flt/dbl, [1] = flt/dbl, [2] = flt/dbl, ...
	               .
	               .
	               .

	 `out_len`
	       [0] -> sigmas1_len[TS_MAX_LEN]
	       [1] -> sigmas2_len[TS_MAX_LEN]

	       `sigmasX_len` (single values)
	             [0] = int, [1] = int, [2] = int, ...  */

	T   (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = (*out)[0];
	T   (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = (*out)[1];
	int (&sigmas1_len)[TS_MAX_LEN]         = (*out_len)[0];
	T pose_out[10];
	for (int i = 0; i < ts_len; i++) {
		sigmas1_len[i] = 0; 

		fn_t(ts[i], pose_out);

		//T &fvalue = pose_out[0]; // double-checked: not used in matlab
		const T &A = pose_out[0], &B = pose_out[1], &C = pose_out[2], 
            &E = pose_out[3], &F = pose_out[4], &G = pose_out[5], &H = pose_out[6],
            &J = pose_out[7], &K = pose_out[8], &L = pose_out[9];

		std::complex<T> delta = sqrt(B*B - 4*A*C);
		std::complex<T> sigma1_m = (-B - delta)/(2*A);
		std::complex<T> sigma1_p = (-B + delta)/(2*A);

		delta = sqrt(F*F - 4*E*G);
		std::complex<T> sigma2_m = (-F - delta)/(2*E);
		std::complex<T> sigma2_p = (-F + delta)/(2*E);

		//% handle case of negative delta
		if (std::abs(std::imag(sigma1_m)) < 1e-4) {
			sigma1_m = std::real(sigma1_m);
			sigma1_p = std::real(sigma1_p);
		} else
			std::cerr << "Ignoring t = " << ts[i] << std::endl;

		if (std::abs(std::imag(sigma2_m)) < 1e-4) {
			sigma2_m = std::real(sigma2_m);
			sigma2_p = std::real(sigma2_p);
		} else
			std::cerr << "Ignoring t = " << ts[i] << std::endl;

		//% Now check to see which pair pass. Only a single pair should pass, in theory.
		//% If not, issue a warning.
		constexpr T my_eps = 1.0;

		if (std::abs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps) {
			sigmas1[i][sigmas1_len[i]++] = sigma1_m.real();
			sigmas2[i][sigmas1_len[i]] = sigma2_m.real();
		}
		if (std::abs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps) {
			if (sigmas1_len[i] != 0) // !isempty(sigmas1[i])
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]++] = sigma1_p.real();
			sigmas2[i][sigmas1_len[i]] = sigma2_m.real();
		}
		if (std::abs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps) {
			if (sigmas1_len[i] != 0)
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]++] = sigma1_p.real();
			sigmas2[i][sigmas1_len[i]] = sigma2_p.real();
		}
		if (std::abs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps) {
			if (sigmas1_len[i] != 0)
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]++] = sigma1_m.real();
			sigmas2[i][sigmas1_len[i]] = sigma2_p.real();
		}
		if (sigmas1_len[i] == 0) // isempty(sigmas1[i])
			std::cerr << "no sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
	}
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
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			output_m[i][j] = 0;
			for (int k = 0; k < 3; k++)
				output_m[i][j] += m1[i][k] * m2[k][j];
		}
}

template<typename T>
void
pose_poly<T>::
get_r_t_from_rhos(
	const int ts_len,
	const T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas1_len)[TS_MAX_LEN],
	const T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN],
	const T (&rhos1)[ROOT_IDS_LEN], const T (&rhos2)[ROOT_IDS_LEN],
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output)[RT_MAX_LEN][4][3], int *output_len
)
{
	T lambdas1[TS_MAX_LEN][TS_MAX_LEN]; T lambdas2[TS_MAX_LEN][TS_MAX_LEN];
	const T DGama[3] = {Gama1[0]-Gama2[0], Gama1[1]-Gama2[1], Gama1[2]-Gama2[2]};
  
	for (int i = 0; i < ts_len; i++) {
    assert(sigmas1_len[i] == sigmas2_len[i]);
    const T dgamas_rhos[3] = {
     rhos1[i]*gama1[0] - rhos2[i]*gama2[0],
     rhos1[i]*gama1[1] - rhos2[i]*gama2[1],
     rhos1[i]*gama1[2] - rhos2[i]*gama2[2]};
		for (int j = 0; j < sigmas1_len[i]; j++) {
			lambdas1[i][j] = 
        (DGama[0]*Tgt1[0]+DGama[1]*Tgt1[1] + DGama[2]*Tgt1[2]) / 
        (dgamas_rhos[0]*(rhos1[i]*tgt1[0] + sigmas1[i][j]*gama1[0]) + 
        dgamas_rhos[1]*(rhos1[i]*tgt1[1] + sigmas1[i][j]*gama1[1]) +
        dgamas_rhos[2]*(rhos1[i]*tgt1[2] + sigmas1[i][j]*gama1[2]));
			lambdas2[i][j] = 
        (DGama[0]*Tgt2[0]+DGama[1]*Tgt2[1] + DGama[2]*Tgt2[2]) /
        (dgamas_rhos[0]*(rhos2[i]*tgt2[0] + sigmas2[i][j]*gama2[0]) + 
        dgamas_rhos[1]*(rhos2[i]*tgt2[1] + sigmas2[i][j]*gama2[1]) +
        dgamas_rhos[2]*(rhos2[i]*tgt2[2] + sigmas2[i][j]*gama2[2]));
		}
	}

	//% Rotation:
	const T A[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2]};
	T inv_A[3][3]; invm3x3(A, inv_A);

	// Matrix containing Rotations and Translations
	T (&RT)[RT_MAX_LEN][4][3] = *output;
	int &RT_len               = *output_len; RT_len = 0;
	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < sigmas1_len[i]; j++, RT_len++) {
			T (&Rots)[4][3] = RT[RT_len]; T (&Transls)[3] = RT[RT_len][3];

			#define B_row(r) \
				rhos1[i]*gama1[(r)] - rhos2[i]*gama2[(r)], \
				lambdas1[i][j]*(rhos1[i]*tgt1[(r)] + sigmas1[i][j]*gama1[(r)]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[(r)] + sigmas2[i][j]*gama2[(r)])

			const T B[3][3] = { B_row(0), B_row(1), B_row(2) };
			multm3x3(B, inv_A, Rots);
      Transls[0] = rhos1[i]*gama1[0] - Rots[0][0] * Gama1[0] - Rots[0][1] * Gama1[1] - Rots[0][2] * Gama1[2];
      Transls[1] = rhos1[i]*gama1[1] - Rots[1][0] * Gama1[0] - Rots[1][1] * Gama1[1] - Rots[1][2] * Gama1[2];
      Transls[2] = rhos1[i]*gama1[2] - Rots[2][0] * Gama1[0] - Rots[2][1] * Gama1[1] - Rots[2][2] * Gama1[2];
		}
	}
}

} // ! namespace P2Pt

#endif // !pose_poly_hxx_


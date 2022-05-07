#ifndef get_sigmas_hxx_
#define get_sigmas_hxx_

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

}

#endif // !get_sigmas_hxx_


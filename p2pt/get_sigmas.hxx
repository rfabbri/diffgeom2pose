#include "pose_from_point_tangents_2_fn_t.hxx"


namespace P2Pt {

template<typename T>
void
pose_poly<T>::
get_sigmas(
	const T (&ts)[ROOT_IDS_LEN], const int ts_len,
	T (*output)[4][SIGMA_LEN][SIGMA_LEN]
)
{
	T my_eps = 1;

	/**************************************************************************
	 `output`
	       [0] -> sigmas1[SIGMA_LEN][SIGMA_LEN]
	       [1] -> sigmas2[SIGMA_LEN][SIGMA_LEN]
	       [2] -> sigmas1_end[SIGMA_LEN][SIGMA_LEN]
	       [3] -> sigmas2_end[SIGMA_LEN][SIGMA_LEN]

	       `sigmasX` (can contain single values or array of values)
	             [0][0] -> float/double
	             [1][0] -> float/double
	             [2][0] -> float/double
	             [3][ ] -> [0] = flt/dbl, [1] = flt/dbl, [2] = flt/dbl, ...
	               .
	               .
	               .

	       `sigmasX_end` (single values)
	             [0][0] -> float/double (cast to integer on use)
	             [0][1] -> float/double (cast to integer on use)
	             [0][2] -> float/double (cast to integer on use)
	               .
	               .
	               .
	 **************************************************************************/

	T (&sigmas1)[SIGMA_LEN][SIGMA_LEN] = (*output)[0];
	T (&sigmas2)[SIGMA_LEN][SIGMA_LEN] = (*output)[1];

	// Stores the number of elements for each row of `sigmasX`
	// Ugly hack, could be set as a struct
	T (&sigmas1_end)[SIGMA_LEN] = (*output)[2][0];
	T (&sigmas2_end)[SIGMA_LEN] = (*output)[3][0];

	static T pose_output[11];

	for (int i = 0; i < ts_len; i++) {
		sigmas1_end[i] = 0;
		sigmas2_end[i] = 0;

		pose_from_point_tangents_2_fn_t(ts[i], &pose_output);

		//T &fvalue = pose_output[0]; // double-checked: not used in matlab
		T &A      = pose_output[1];
		T &B      = pose_output[2];
		T &C      = pose_output[3];
		T &E      = pose_output[4];
		T &F      = pose_output[5];
		T &G      = pose_output[6];
		T &H      = pose_output[7];
		T &J      = pose_output[8];
		T &K      = pose_output[9];
		T &L      = pose_output[10];

		// TODO: Optimize quadratic function algorithm
		std::complex<T> delta1 = sqrt(B*B - 4*A*C);
		std::complex<T> sigma1_m = (-B - delta1)/(2*A);
		std::complex<T> sigma1_p = (-B + delta1)/(2*A);

		std::complex<T> delta2 = sqrt(F*F - 4*E*G);
		std::complex<T> sigma2_m = (-F - delta2)/(2*E);
		std::complex<T> sigma2_p = (-F + delta2)/(2*E);

		//% handle case of negative delta
		if (std::abs(std::imag(sigma1_m)) < 1e-4) {
			sigma1_m = std::real(sigma1_m);
			sigma1_p = std::real(sigma1_p);
		} else {
			std::cerr << "Ignoring t = " << ts[i] << std::endl;
		}

		if (std::abs(std::imag(sigma2_m)) < 1e-4) {
			sigma2_m = std::real(sigma2_m);
			sigma2_p = std::real(sigma2_p);
		} else {
			std::cerr << "Ignoring t = " << ts[i] << std::endl;
		}

		//% Now check to see which pair pass. Only a single pair should pass, in theory.
		//% If not, issue a warning.

		if (std::abs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps) {
			sigmas1[i][(int)sigmas1_end[i]++] = sigma1_m.real();
			sigmas2[i][(int)sigmas2_end[i]++] = sigma2_m.real();
		}

		if (std::abs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps) {
			if (sigmas1_end[i] != 0) // !isempty(sigmas1[i])
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][(int)sigmas1_end[i]++] = sigma1_p.real();
			sigmas2[i][(int)sigmas2_end[i]++] = sigma2_m.real();
		}

		if (std::abs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps) {
			if (sigmas1_end[i] != 0)
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][(int)sigmas1_end[i]++] = sigma1_p.real();
			sigmas2[i][(int)sigmas2_end[i]++] = sigma2_p.real();
		}

		if (std::abs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps) {
			if (sigmas1_end[i] != 0)
				std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][(int)sigmas1_end[i]++] = sigma1_m.real();
			sigmas2[i][(int)sigmas2_end[i]++] = sigma2_p.real();
		}

		if (sigmas1_end[i] == 0) // isempty(sigmas1[i])
			std::cerr << "no sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;

	}
}

}


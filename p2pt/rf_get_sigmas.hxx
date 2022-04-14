namespace P2Pt {

template<typename T>
void
rf_get_sigmas(T ts[t_vector_len], static int ts_len, T output[11]) // TODO: Fix size of output
{
	T my_eps = 1;

	// TODO: The size of this array can possibly be reduced.
	// Figure out the maximum number of 1s in `root_ids[]`
	static T sigmas1[ts_len][ts_len];
	static T sigmas2[ts_len][ts_len];

	// Stores the number of elements for each row of `sigmasX`
	static int end_sigmas1[ts_len];
	static int end_sigmas2[ts_len];

	static T rf_pose_output[11];

	for (int i = 0; i < ts_len; i++) {
		end_sigmas1[i] = 0;
		end_sigmas2[i] = 0;

		rf_pose_from_point_tangents_2_fn_t(ts[i], rf_pose_output);

		// Getting value instead of pointer to avoid confusion between
		// dereference `*` and multiplication `*` in the formulas below
		static T fvalue = rf_pose_output[0];
		static T A      = rf_pose_output[1];
		static T B      = rf_pose_output[2];
		static T C      = rf_pose_output[3];
		static T E      = rf_pose_output[4];
		static T F      = rf_pose_output[5];
		static T G      = rf_pose_output[6];
		static T H      = rf_pose_output[7];
		static T J      = rf_pose_output[8];
		static T K      = rf_pose_output[9];
		static T L      = rf_pose_output[10];

		// TODO: Optimize quadratic function algorithm
		// TODO: Deal with complex results
		// TODO: use std::complex; maybe use another approach later
		T delta1 = sqrt(B*B - 4*A*C);
		T sigma1_m = (-B - delta1)/(2*A);
		T sigma2_p = (-B + delta1)/(2*A);

		T delta2 = sqrt(F*F - 4*E*G);
		T sigma1_m = (-F - delta2)/(2*E);
		T sigma2_p = (-F + delta22/(2*E);

		//% handle case of negative delta
		// TODO: Get real/imaginary part of answer using some std library function
		// TODO: Check if this precision can be turned into a constant
		if (abs(imag(sigma1_m)) < 1e-4) {
			sigma1_m = real(sigma1_m);
			sigma1_p = real(sigma1_p);
		} else {
			// TODO: Check if this needs to be redirected to stderr
			std::cout << "Ignoring t = " << ts[i] << std::endl;
		}

		if (abs(imag(sigma2_m)) < 1e-4) {
			sigma2_m = real(sigma2_m);
			sigma2_p = real(sigma2_p);
		} else {
			std::cout << "Ignoring t = " << ts[i] << std::endl;
		}

		//% Now check to see which pair pass. Only a single pair should pass, in theory.
		//% If not, issue a warning.

		if (abs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps) {
			sigmas1[i][end_sigmas1[i]++] = sigma1_m;
			sigmas2[i][end_sigmas2[i]++] = sigma2_m;
		}

		if (abs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps) {
			// TODO: Mark `sigmas1` and `sigmas2` indexes as empty or not. Nullptr? Or check end1 == 0
			if (end_sigmas1 != 0) // !isempty(sigmas1[i])
				std::cout << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][end_sigmas1[i]++] = sigma1_p;
			sigmas2[i][end_sigmas2[i]++] = sigma2_m;
		}

		if (abs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps) {
			if (end_sigmas1 != 0)
				std::cout << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][end_sigmas1[i]++] = sigma1_p;
			sigmas2[i][end_sigmas2[i]++] = sigma2_p;
		}

		if (abs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps) {
			if (end_sigmas1 != 0)
				std::cout << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][end_sigmas1[i]++] = sigma1_m;
			sigmas2[i][end_sigmas2[i]++] = sigma2_p;
		}

		if (end_sigmas1 == 0) // isempty(sigmas1[i])
			std::cout << "no sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;

	}
}

}


#ifndef poly_h_
#define poly_h_

constexpr int T_VECTOR_LEN = 2001;
constexpr int ROOT_IDS_LEN = T_VECTOR_LEN - 1;

// TODO: THIS IS TEMPORARY! FOR TEST ONLY
constexpr int TS_LEN = 4;
constexpr int RT_LEN = 7;

namespace P2Pt {

template<typename T>
struct pose_poly {
	T A0, A1, A2;
	T B0, B1, B2, B3;
	T C0, C1, C2, C3, C4;
	T E0, E1, E2;
	T F0, F1, F2, F3;
	T G0, G1, G2, G3, G4;
	T H0, H1, H2, H3, H4;
	T J0, J1, J2, J3;
	T K0, K1, K2, K3;
	T L0, L1, L2;
	T alpha, beta, theta;
	// TODO: Add cos(theta) and sin(theta) as members
	void pose_from_point_tangents_2(
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3]
	);
	void find_bounded_root_intervals(
		const T (&t_vector)[T_VECTOR_LEN], T (*root_ids_output)[ROOT_IDS_LEN]
	);
	void sample_pose_poly(
		const T (&t)[T_VECTOR_LEN], T (*output)[11][T_VECTOR_LEN] // = {fvalue, A, B, C, E, F, G, H, J, K, L}
	);
	void pose_from_point_tangents_2_fn_t_for_root(
		const T t, T (*output)[11] // = {fvalue, A, B, C, E, F, G, H, J, K, L}
	);
	void pose_from_point_tangents_2_fn_t(
		const T t, T (*output)[11] // = {fvalue, A, B, C, E, F, G, H, J, K, L}
	);

	/**************************************************************************
	 THINGS THAT NEED TO BE REFACTORED (using len, row/col arguments)
	 - `rhos_from_root_ids`:
	     + `output` [OK]
	 - `get_sigmas`:
	     + `ts`     [OK]
	     + `output` [OK]
	 - `get_r_t_from_rhos`
	     + `sigmas1`
		 + `sigmas2`
		 + `rhos1`
		 + `rhos2`
	 **************************************************************************/

	void rhos_from_root_ids(
		const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
		T (*output)[8][ROOT_IDS_LEN] // = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts, ts_len}
	);
	void get_sigmas(
		const int ts_len,
		const T ts[],
		T *output // output[4][sigma_len][sigma_len] = {sigmas1, sigmas2, sigmas1_end, sigmas2_end}
	);
	void get_r_t_from_rhos(
		const int ts_len,
		// (&sigmas1)[ts_len][ts_len]
		// (&sigmas2)[ts_len][ts_len]
		// (&sigmas1_end)[ts_len]
		// (&sigmas2_end)[ts_len]
		// (&rhos1)[ts_len]
		// (&rhos2)[ts_len]
		// (*output)[RT_LEN][4][3]
		const T sigmas1[], const T sigmas1_end[],
		const T sigmas2[], const T sigmas2_end[],
		const T rhos1[], const T rhos2[],
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3],
		T *output
	);
};

}

#endif // poly_h_

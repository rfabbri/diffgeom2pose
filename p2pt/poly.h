#ifndef poly_h_
#define poly_h_

constexpr int T_VECTOR_LEN = 2001;
constexpr int ROOT_IDS_LEN = T_VECTOR_LEN - 1;

// TODO: Assume a reasonable length for `ts`. Check if it can be longer than 4.
constexpr int TS_MAX_LEN = 4;
//constexpr int SIGMAS_MAX_LEN = (TS_MAX_LEN * TS_MAX_LEN);
constexpr int RT_MAX_LEN = (TS_MAX_LEN * TS_MAX_LEN);


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
		const T (&t)[T_VECTOR_LEN], T (*output)[11][T_VECTOR_LEN]
	);
	void pose_from_point_tangents_2_fn_t_for_root(
		const T t, T (*output)[11]
	);
	void pose_from_point_tangents_2_fn_t(
		const T t, T (*output)[11]
	);
	void rhos_from_root_ids(
		const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
		T (*output)[8][ROOT_IDS_LEN] /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts, ts_len} */
	);
	void get_sigmas(
		const int ts_len, const T (&ts)[ROOT_IDS_LEN],
		T (*output)[4][TS_MAX_LEN][TS_MAX_LEN]
	);
	void get_r_t_from_rhos(
		const int ts_len,
		const T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN], const T (&sigmas1_end)[TS_MAX_LEN],
		const T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN], const T (&sigmas2_end)[TS_MAX_LEN],
		const T (&rhos1)[ROOT_IDS_LEN], const T (&rhos2)[ROOT_IDS_LEN],
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3],
		T (*output)[RT_MAX_LEN + 1][4][3]
	);
};

}

#endif // poly_h_

#ifndef poly_h_
#define poly_h_

// TODO: Make this global across files
constexpr int T_VECTOR_LEN = 2001;
constexpr int ROOT_IDS_LEN = T_VECTOR_LEN - 1;

// TODO: THIS IS TEMPORARY! FOR TEST ONLY
constexpr int TS_LEN = 4;

// TODO: Figure max length of sigma. Should be equal to `TS_LEN`?
constexpr int SIGMA_LEN = 10;

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
		const T gama1[3], const T tgt1[3],
		const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3],
		const T Gama2[3], const T Tgt2[3]
	);
	void find_bounded_root_intervals(const T t_vector[T_VECTOR_LEN], T root_ids[ROOT_IDS_LEN]);
	void sample_pose_poly(const T t[T_VECTOR_LEN], T output[11][T_VECTOR_LEN]);
	void pose_from_point_tangents_2_fn_t_for_root(const T t, T output[11]);
	void pose_from_point_tangents_2_fn_t(const T t, T output[11]);
	void rhos_from_root_ids(
		const T t_vector[T_VECTOR_LEN], const T root_ids[T_VECTOR_LEN],
		T output[7][T_VECTOR_LEN] /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts} */
	);
	void get_sigmas(
		const T ts[T_VECTOR_LEN], const int ts_len,
		T output[4][SIGMA_LEN][SIGMA_LEN]
	);
	void get_r_t_from_rhos(
		const int ts_len,
		T sigmas1[TS_LEN][TS_LEN], int end_sigmas1[TS_LEN],
		T sigmas2[TS_LEN][TS_LEN], int end_sigmas2[TS_LEN],
		T rhos1[T_VECTOR_LEN], T rhos2[T_VECTOR_LEN]
	);
};

}

#endif // poly_h_

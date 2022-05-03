#ifndef poly_h_
#define poly_h_

constexpr int T_VECTOR_LEN = 2001;
constexpr int ROOT_IDS_LEN = T_VECTOR_LEN - 1;

// TODO: Assume a reasonable length for `ts`. Check if it can be longer than 8.
constexpr int TS_MAX_LEN = 8;
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

	void pose_from_point_tangents_2(
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3]
	);
	void find_bounded_root_intervals(
		const T (&t_vector)[T_VECTOR_LEN], T (*root_ids_output)[ROOT_IDS_LEN]
	);
	T pose_from_point_tangents_2_fn_t(
		const T t, T (*output)[10] = nullptr
	);
	void rhos_from_root_ids(
		const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
		T (*output)[3][ROOT_IDS_LEN], int *output_ts_len
	);
	void get_sigmas(
		const int ts_len, const T (&ts)[ROOT_IDS_LEN],
		T (*output)[2][TS_MAX_LEN][TS_MAX_LEN], int (*output_len)[2][TS_MAX_LEN]
	);
	void get_r_t_from_rhos(
		const int ts_len,
		const T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas1_len)[TS_MAX_LEN],
		const T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas2_len)[TS_MAX_LEN],
		const T (&rhos1)[ROOT_IDS_LEN], const T (&rhos2)[ROOT_IDS_LEN],
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3],
		T (*output)[RT_MAX_LEN][4][3], int *output_len
	);
	T operator()(T t) {
		return pose_poly<T>::pose_from_point_tangents_2_fn_t(t);
	}
};

}

#endif // !poly_h_


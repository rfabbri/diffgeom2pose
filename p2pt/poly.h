#ifndef poly_h_
#define poly_h_

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
	void rf_pose_from_point_tangents_2(
		const T gama1[3], const T tgt1[3],
		const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3],
		const T Gama2[3], const T Tgt2[3]
	);
	void rf_find_bounded_root_intervals(const T t_vector[2001], T root_ids[2001]);
	void rf_sample_pose_poly(
		const T t[2001],
		T A[2001], T B[2001], T C[2001], T E[2001], T F[2001],
		T G[2001], T H[2001], T J[2001], T K[2001], T L[2001], T fvalue[2001]
	);
};

}

#endif // poly_h_

#ifndef pose_poly_h_
#define pose_poly_h_

namespace P2Pt {

constexpr int T_VECTOR_LEN = 2001;
constexpr int ROOT_IDS_LEN = T_VECTOR_LEN - 1;

// TODO: Assume a reasonable length for `ts`. Check if it can be longer than 8.
constexpr int TS_MAX_LEN = 8;
constexpr int RT_MAX_LEN = (TS_MAX_LEN * TS_MAX_LEN);

template<typename T>
struct pose_poly {
	T A0, A1, A2,
		B0, B1, B2, B3,
		C0, C1, C2, C3, C4,
		E0, E1, E2,
		F0, F1, F2, F3,
		G0, G1, G2, G3, G4,
		H0, H1, H2, H3, H4,
		J0, J1, J2, J3,
		K0, K1, K2, K3,
		L0, L1, L2,
		alpha, beta, theta,
    sth, cth;

	void pose_from_point_tangents_2(
		const T (&gama1)[3], const T (&tgt1)[3],
		const T (&gama2)[3], const T (&tgt2)[3],
		const T (&Gama1)[3], const T (&Tgt1)[3],
		const T (&Gama2)[3], const T (&Tgt2)[3]
	);
  
	inline void find_bounded_root_intervals(
		const T (&t_vector)[T_VECTOR_LEN], T (*root_ids_output)[ROOT_IDS_LEN]
	)
  {
    T curr_val = fn_t(t_vector[0]), next_val;
    for (int i = 0; i < ROOT_IDS_LEN; i++) {
      next_val = fn_t(t_vector[i+1]);
      (*root_ids_output)[i] = (curr_val * next_val) < 0;
      curr_val = next_val;
    }
  }
  
  // o: output
	inline T fn_t(const T t, T o[10]) 
  {
    T &A = o[0]; T &B = o[1]; T &C = o[2]; T &E = o[3]; 
    T &F = o[4]; T &G = o[5]; T &H = o[6]; T &J = o[7]; 
    T &K = o[8]; T &L = o[9];

    //%function of t part :
    const T t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t, t6 = t5*t, t7 = t6*t, t8 = t7*t;

    // `(t^2 + 1)` integer powers
    // TODO(check) when removing static, the tests accuse loss of precision
    const T t2p12 = (t2 + 1.) * (t2 + 1.), t2p13 = t2p12 * (t2 + 1.), t2p14 = t2p13 * (t2 + 1.);

    A = (A0+A1*t+A2*t2-A1*t3+A0*t4)/t2p12;
    B = (B0+B1*t+B2*t2+B3*t3-B2*t4+B1*t5-B0*t6)/t2p13;
    C = (C0+C1*t+C2*t2+C3*t3+C4*t4-C3*t5+C2*t6-C1*t7+C0*t8)/t2p14;
    E = (E0+E1*t+E2*t2-E1*t3+E0*t4)/t2p12;
    F = (F0+F1*t+F2*t2+F3*t3-F2*t4+F1*t5-F0*t6)/t2p13;
    G = (G0+G1*t+G2*t2+G3*t3+G4*t4-G3*t5+G2*t6-G1*t7+G0*t8)/t2p14;
    H = (H0+H1*t+H2*t2+H3*t3+H4*t4-H3*t5+H2*t6-H1*t7+H0*t8)/t2p14;
    J = (J0+J1*t+J2*t2+J3*t3-J2*t4+J1*t5-J0*t6)/t2p13;
    K = (K0+K1*t+K2*t2+K3*t3-K2*t4+K1*t5-K0*t6)/t2p13;
    L = (L0+L1*t+L2*t2-L1*t3+L0*t4)/t2p12;

    const T
    A2=A*A, B2=B*B, C2=C*C, E2=E*E, F2=F*F, G2=G*G, H2=H*H, 
    H3=H2*H, H4=H3*H, J2=J*J, J3=J2*J, K2=K*K, K3=K2*K, L2=L*L, L3=L2*L;

    return // fvalue_terms[X]
    // TODO: Analyze the need of extra/lower precison for the calculation of these terms.
    E2*B2*H2*J2 +G2*C2*L*L*L*L +G2*A2*K3*K +E2*A2*H4 +E2*C2*J3*J
    +-2.*E*A*H2*G*C*L2 +2.*E2*A*H2*C*J2 +-2.*E2*C*J3*B*H +2.*E*C2*J2*G*L2
    +2.*E*A2*H2*G*K2 +-2.*E2*A*H3*B*J +-2.*E*A*H2*G*B*K*L +-2.*E*C*J2*G*B*K*L
    +-2.*E*C*J2*G*A*K2 +-2.*E*B*H*J*G*C*L2 +-2.*E*B*H*J*G*A*K2 +G2*B2*K2*L2
    +-2.*G2*B*K*L3*C +-2.*G2*B*K3*L*A +2.*G2*C*L2*A*K2 +-2.*F*E*A2*H3*K
    +-2.*F*E*A*H*K*C*J2 +3.*F*E*A*H2*K*B*J +3.*F*A*H*K2*G*B*L
    +-2.*F*A*H*K*G*C*L2 +-2.*F*A2*H*K3*G +F*E*B*H3*L*A +3.*F*E*B*H*L*C*J2
    +-1.*F*E*B2*H2*L*J +-1.*F*B2*H*L2*G*K +F*B*H*L3*G*C +F*E*B*K*J3*C
    +-1.*F*E*B2*K*J2*H +-1.*F*B2*K2*J*G*L +3.*F*B*K*J*G*C*L2 +F*B*K3*J*G*A
    +-2.*F*E*C*J*L*A*H2 +-2.*F*E*C2*J3*L +-2.*F*C2*J*L3*G +-2.*F*C*J*L*G*A*K2
    +F2*A2*K2*H2 +F2*A*K2*C*J2 +-1.*F2*A*K2*B*H*J +-1.*F2*B*K*L*A*H2
    +-1.*F2*B*K*L*C*J2 +F2*B2*K*L*H*J +F2*C*L2*A*H2 +F2*C2*L2*J2
    +-1.*F2*C*L2*B*H*J +G*E*B2*H2*L2 +G*E*B2*K2*J2 +8.*G*E*A*H*K*C*J*L;
  }
  
	inline T fn_t(const T t) { T b[10]; return fn_t(t, b); }
	inline T operator()(T t) { return fn_t(t); }
  
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
};

}

#endif // !poly_h_


#ifndef pose_from_point_tangents_2_fn_t_hxx_
#define pose_from_point_tangents_2_fn_t_hxx_

namespace P2Pt {

template<typename T>
inline T
pose_poly<T>::
pose_from_point_tangents_2_fn_t(const T t, T (*output)[10] /* = nullptr */)
{
	static T buf[10];

	//%function of t part :
	T& A = *output ? (*output)[0] : buf[0];
	T& B = *output ? (*output)[1] : buf[1];
	T& C = *output ? (*output)[2] : buf[2];
	T& E = *output ? (*output)[3] : buf[3];
	T& F = *output ? (*output)[4] : buf[4];
	T& G = *output ? (*output)[5] : buf[5];
	T& H = *output ? (*output)[6] : buf[6];
	T& J = *output ? (*output)[7] : buf[7];
	T& K = *output ? (*output)[8] : buf[8];
	T& L = *output ? (*output)[9] : buf[9];

	// `t` integer powers
	static T t_pow2, t_pow3, t_pow4, t_pow5, t_pow6, t_pow7, t_pow8;

	t_pow2 = t * t;
	t_pow3 = t_pow2 * t;
	t_pow4 = t_pow3 * t;
	t_pow5 = t_pow4 * t;
	t_pow6 = t_pow5 * t;
	t_pow7 = t_pow6 * t;
	t_pow8 = t_pow7 * t;

	// `(t^2 + 1)` integer powers
	static T t_pow2_plus1_pow2, t_pow2_plus1_pow3, t_pow2_plus1_pow4;

	t_pow2_plus1_pow2 = (t_pow2 + 1) * (t_pow2 + 1);
	t_pow2_plus1_pow3 = t_pow2_plus1_pow2 * (t_pow2 + 1);
	t_pow2_plus1_pow4 = t_pow2_plus1_pow3 * (t_pow2 + 1);

	// Denominators
	T &A_den = t_pow2_plus1_pow2;
	T &B_den = t_pow2_plus1_pow3;
	T &C_den = t_pow2_plus1_pow4;
	T &E_den = t_pow2_plus1_pow2;
	T &F_den = t_pow2_plus1_pow3;
	T &G_den = t_pow2_plus1_pow4;
	T &H_den = t_pow2_plus1_pow4;
	T &J_den = t_pow2_plus1_pow3;
	T &K_den = t_pow2_plus1_pow3;
	T &L_den = t_pow2_plus1_pow2;

	// Calculations
	A = (A0 + A1 * t + A2 * t_pow2 - A1 * t_pow3 + A0 * t_pow4) / A_den;
	B = (B0 + B1 * t + B2 * t_pow2 + B3 * t_pow3 - B2 * t_pow4 + B1 * t_pow5 - B0 * t_pow6) / B_den;
	C = (C0 + C1 * t + C2 * t_pow2 + C3 * t_pow3 + C4 * t_pow4 - C3 * t_pow5 + C2 * t_pow6 - C1 * t_pow7 + C0 * t_pow8) / C_den;
	E = (E0 + E1 * t + E2 * t_pow2 - E1 * t_pow3 + E0 * t_pow4) / E_den;
	F = (F0 + F1 * t + F2 * t_pow2 + F3 * t_pow3 - F2 * t_pow4 + F1 * t_pow5 - F0 * t_pow6) / F_den;
	G = (G0 + G1 * t + G2 * t_pow2 + G3 * t_pow3 + G4 * t_pow4 - G3 * t_pow5 + G2 * t_pow6 - G1 * t_pow7 + G0 * t_pow8) / G_den;
	H = (H0 + H1 * t + H2 * t_pow2 + H3 * t_pow3 + H4 * t_pow4 - H3 * t_pow5 + H2 * t_pow6 - H1 * t_pow7 + H0 * t_pow8) / H_den;
	J = (J0 + J1 * t + J2 * t_pow2 + J3 * t_pow3 - J2 * t_pow4 + J1 * t_pow5 - J0 * t_pow6) / J_den;
	K = (K0 + K1 * t + K2 * t_pow2 + K3 * t_pow3 - K2 * t_pow4 + K1 * t_pow5 - K0 * t_pow6) / K_den;
	L = (L0 + L1 * t + L2 * t_pow2 - L1 * t_pow3 + L0 * t_pow4) / L_den;

	static T A_pow2, B_pow2, C_pow2, E_pow2, F_pow2, G_pow2, H_pow2, J_pow2, K_pow2, L_pow2;
	static T H_pow3, J_pow3, K_pow3, L_pow3;
	static T H_pow4, J_pow4, K_pow4, L_pow4;

	// 2nd power          3rd power               4th power
	A_pow2 = A * A;
	B_pow2 = B * B;
	C_pow2 = C * C;
	E_pow2 = E * E;
	F_pow2 = F * F;
	G_pow2 = G * G;
	H_pow2 = H * H;    H_pow3 = H_pow2 * H;    H_pow4 = H_pow3 * H;
	J_pow2 = J * J;    J_pow3 = J_pow2 * J;    J_pow4 = J_pow3 * J;
	K_pow2 = K * K;    K_pow3 = K_pow2 * K;    K_pow4 = K_pow3 * K;
	L_pow2 = L * L;    L_pow3 = L_pow2 * L;    L_pow4 = L_pow3 * L;

	static T fvalue;
	fvalue = 0;

	// TODO: Analyze the need of extra/lower precison for the calculation of these terms.
	// Report about it in the thesis
	// fvalue_terms[X]
	/*  0 */ fvalue +=      E_pow2 * B_pow2 * H_pow2 * J_pow2;
	/*  1 */ fvalue +=      G_pow2 * C_pow2 * L_pow4;
	/*  2 */ fvalue +=      G_pow2 * A_pow2 * K_pow4;
	/*  3 */ fvalue +=      E_pow2 * A_pow2 * H_pow4;
	/*  4 */ fvalue +=      E_pow2 * C_pow2 * J_pow4;
	/*  5 */ fvalue += -2 * E      * A      * H_pow2 * G      * C      * L_pow2;
	/*  6 */ fvalue +=  2 * E_pow2 * A      * H_pow2 * C      * J_pow2;
	/*  7 */ fvalue += -2 * E_pow2 * C      * J_pow3 * B      * H;
	/*  8 */ fvalue +=  2 * E      * C_pow2 * J_pow2 * G      * L_pow2;
	/*  9 */ fvalue +=  2 * E      * A_pow2 * H_pow2 * G      * K_pow2;
	/* 10 */ fvalue += -2 * E_pow2 * A      * H_pow3 * B      * J;
	/* 11 */ fvalue += -2 * E      * A      * H_pow2 * G      * B      * K      * L;
	/* 12 */ fvalue += -2 * E      * C      * J_pow2 * G      * B      * K      * L;
	/* 13 */ fvalue += -2 * E      * C      * J_pow2 * G      * A      * K_pow2;
	/* 14 */ fvalue += -2 * E      * B      * H      * J      * G      * C      * L_pow2;
	/* 15 */ fvalue += -2 * E      * B      * H      * J      * G      * A      * K_pow2;
	/* 16 */ fvalue +=      G_pow2 * B_pow2 * K_pow2 * L_pow2;
	/* 17 */ fvalue += -2 * G_pow2 * B      * K      * L_pow3 * C;
	/* 18 */ fvalue += -2 * G_pow2 * B      * K_pow3 * L      * A;
	/* 19 */ fvalue +=  2 * G_pow2 * C      * L_pow2 * A      * K_pow2;
	/* 20 */ fvalue += -2 * F      * E      * A_pow2 * H_pow3 * K;
	/* 21 */ fvalue += -2 * F      * E      * A      * H      * K      * C      * J_pow2;
	/* 22 */ fvalue +=  3 * F      * E      * A      * H_pow2 * K      * B      * J;
	/* 23 */ fvalue +=  3 * F      * A      * H      * K_pow2 * G      * B      * L;
	/* 24 */ fvalue += -2 * F      * A      * H      * K      * G      * C      * L_pow2;
	/* 25 */ fvalue += -2 * F      * A_pow2 * H      * K_pow3 * G;
	/* 26 */ fvalue +=      F      * E      * B      * H_pow3 * L      * A;
	/* 27 */ fvalue +=  3 * F      * E      * B      * H      * L      * C      * J_pow2;
	/* 28 */ fvalue += -1 * F      * E      * B_pow2 * H_pow2 * L      * J;
	/* 29 */ fvalue += -1 * F      * B_pow2 * H      * L_pow2 * G      * K;
	/* 30 */ fvalue +=      F      * B      * H      * L_pow3 * G      * C;
	/* 31 */ fvalue +=      F      * E      * B      * K      * J_pow3 * C;
	/* 32 */ fvalue += -1 * F      * E      * B_pow2 * K      * J_pow2 * H;
	/* 33 */ fvalue += -1 * F      * B_pow2 * K_pow2 * J      * G      * L;
	/* 34 */ fvalue +=  3 * F      * B      * K      * J      * G      * C      * L_pow2;
	/* 35 */ fvalue +=      F      * B      * K_pow3 * J      * G      * A;
	/* 36 */ fvalue += -2 * F      * E      * C      * J      * L      * A      * H_pow2;
	/* 37 */ fvalue += -2 * F      * E      * C_pow2 * J_pow3 * L;
	/* 38 */ fvalue += -2 * F      * C_pow2 * J      * L_pow3 * G;
	/* 39 */ fvalue += -2 * F      * C      * J      * L      * G      * A      * K_pow2;
	/* 40 */ fvalue +=      F_pow2 * A_pow2 * K_pow2 * H_pow2;
	/* 41 */ fvalue +=      F_pow2 * A      * K_pow2 * C      * J_pow2;
	/* 42 */ fvalue += -1 * F_pow2 * A      * K_pow2 * B      * H      * J;
	/* 43 */ fvalue += -1 * F_pow2 * B      * K      * L      * A      * H_pow2;
	/* 44 */ fvalue += -1 * F_pow2 * B      * K      * L      * C      * J_pow2;
	/* 45 */ fvalue +=      F_pow2 * B_pow2 * K      * L      * H      * J;
	/* 46 */ fvalue +=      F_pow2 * C      * L_pow2 * A      * H_pow2;
	/* 47 */ fvalue +=      F_pow2 * C_pow2 * L_pow2 * J_pow2;
	/* 48 */ fvalue += -1 * F_pow2 * C      * L_pow2 * B      * H      * J;
	/* 49 */ fvalue +=      G      * E      * B_pow2 * H_pow2 * L_pow2;
	/* 50 */ fvalue +=      G      * E      * B_pow2 * K_pow2 * J_pow2;
	/* 51 */ fvalue +=  8 * G      * E      * A      * H      * K      * C      * J      * L;

	return fvalue;
}

}

#endif // !pose_from_point_tangents_2_fn_t_hxx_


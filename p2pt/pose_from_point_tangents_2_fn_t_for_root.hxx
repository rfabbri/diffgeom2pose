#include "poly.h"

namespace P2Pt {

template<typename T>
T
pose_poly<T>::
pose_from_point_tangents_2_fn_t_for_root(const T t)
{
	// TODO: See if this can be reused from `sample_pose_poly`, even though this function
	// only applies to scalar values rather than vectors. Function overloading?

	//% perturb = -4.9e6;
	//% perturb = -1.7e8;
	//% perturb = -2.55e7;
	//static T perturb;
	//perturb = 0;
	//% same as  pose_from_point_tangents_2_fn_t, but polynomial value might be
	//% perturbed by max of(perturb) so that roots appear for noisy data.This wasn't
	//% used in the end.
	//%

	//%function of t part :
	static T fvalue, A, B, C, E, F, G, H, J, K, L;

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

	// `output`
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

	static T fvalue_terms[52];

	fvalue_terms[0]  =      E_pow2 * B_pow2 * H_pow2 * J_pow2;
	fvalue_terms[1]  =      G_pow2 * C_pow2 * L_pow4;
	fvalue_terms[2]  =      G_pow2 * A_pow2 * K_pow4;
	fvalue_terms[3]  =      E_pow2 * A_pow2 * H_pow4;
	fvalue_terms[4]  =      E_pow2 * C_pow2 * J_pow4;
	fvalue_terms[5]  = -2 * E      * A      * H_pow2 * G      * C      * L_pow2;
	fvalue_terms[6]  =  2 * E_pow2 * A      * H_pow2 * C      * J_pow2;
	fvalue_terms[7]  = -2 * E_pow2 * C      * J_pow3 * B      * H;
	fvalue_terms[8]  =  2 * E      * C_pow2 * J_pow2 * G      * L_pow2;
	fvalue_terms[9]  =  2 * E      * A_pow2 * H_pow2 * G      * K_pow2;
	fvalue_terms[10] = -2 * E_pow2 * A      * H_pow3 * B      * J;
	fvalue_terms[11] = -2 * E      * A      * H_pow2 * G      * B      * K      * L;
	fvalue_terms[12] = -2 * E      * C      * J_pow2 * G      * B      * K      * L;
	fvalue_terms[13] = -2 * E      * C      * J_pow2 * G      * A      * K_pow2;
	fvalue_terms[14] = -2 * E      * B      * H      * J      * G      * C      * L_pow2;
	fvalue_terms[15] = -2 * E      * B      * H      * J      * G      * A      * K_pow2;
	fvalue_terms[16] =      G_pow2 * B_pow2 * K_pow2 * L_pow2;
	fvalue_terms[17] = -2 * G_pow2 * B      * K      * L_pow3 * C;
	fvalue_terms[18] = -2 * G_pow2 * B      * K_pow3 * L      * A;
	fvalue_terms[19] =  2 * G_pow2 * C      * L_pow2 * A      * K_pow2;
	fvalue_terms[20] = -2 * F      * E      * A_pow2 * H_pow3 * K;
	fvalue_terms[21] = -2 * F      * E      * A      * H      * K      * C      * J_pow2;
	fvalue_terms[22] =  3 * F      * E      * A      * H_pow2 * K      * B      * J;
	fvalue_terms[23] =  3 * F      * A      * H      * K_pow2 * G      * B      * L;
	fvalue_terms[24] = -2 * F      * A      * H      * K      * G      * C      * L_pow2;
	fvalue_terms[25] = -2 * F      * A_pow2 * H      * K_pow3 * G;
	fvalue_terms[26] =      F      * E      * B      * H_pow3 * L      * A;
	fvalue_terms[27] =  3 * F      * E      * B      * H      * L      * C      * J_pow2;
	fvalue_terms[28] = -1 * F      * E      * B_pow2 * H_pow2 * L      * J;
	fvalue_terms[29] = -1 * F      * B_pow2 * H      * L_pow2 * G      * K;
	fvalue_terms[30] =      F      * B      * H      * L_pow3 * G      * C;
	fvalue_terms[31] =      F      * E      * B      * K      * J_pow3 * C;
	fvalue_terms[32] = -1 * F      * E      * B_pow2 * K      * J_pow2 * H;
	fvalue_terms[33] = -1 * F      * B_pow2 * K_pow2 * J      * G      * L;
	fvalue_terms[34] =  3 * F      * B      * K      * J      * G      * C      * L_pow2;
	fvalue_terms[35] =      F      * B      * K_pow3 * J      * G      * A;
	fvalue_terms[36] = -2 * F      * E      * C      * J      * L      * A      * H_pow2;
	fvalue_terms[37] = -2 * F      * E      * C_pow2 * J_pow3 * L;
	fvalue_terms[38] = -2 * F      * C_pow2 * J      * L_pow3 * G;
	fvalue_terms[39] = -2 * F      * C      * J      * L      * G      * A      * K_pow2;
	fvalue_terms[40] =      F_pow2 * A_pow2 * K_pow2 * H_pow2;
	fvalue_terms[41] =      F_pow2 * A      * K_pow2 * C      * J_pow2;
	fvalue_terms[42] = -1 * F_pow2 * A      * K_pow2 * B      * H      * J;
	fvalue_terms[43] = -1 * F_pow2 * B      * K      * L      * A      * H_pow2;
	fvalue_terms[44] = -1 * F_pow2 * B      * K      * L      * C      * J_pow2;
	fvalue_terms[45] =      F_pow2 * B_pow2 * K      * L      * H      * J;
	fvalue_terms[46] =      F_pow2 * C      * L_pow2 * A      * H_pow2;
	fvalue_terms[47] =      F_pow2 * C_pow2 * L_pow2 * J_pow2;
	fvalue_terms[48] = -1 * F_pow2 * C      * L_pow2 * B      * H      * J;
	fvalue_terms[49] =      G      * E      * B_pow2 * H_pow2 * L_pow2;
	fvalue_terms[50] =      G      * E      * B_pow2 * K_pow2 * J_pow2;
	fvalue_terms[51] =  8 * G      * E      * A      * H      * K      * C      * J      * L;

	fvalue = 0;
	for (int i = 0; i < 52; i++)
		fvalue += fvalue_terms[i] /* + perturb */;
	return fvalue;
}

}


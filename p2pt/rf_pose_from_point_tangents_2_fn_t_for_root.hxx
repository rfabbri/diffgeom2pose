#include "poly.h"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
rf_pose_from_point_tangents_2_fn_t_for_root(const T t, T output[11])
{
	// TODO: See if this can be reused from `rf_sample_pose_poly`, even though this function
	// only applies to scalar values rather than vectors. Function overloading?

	//% perturb = -4.9e6;
	//% perturb = -1.7e8;
	//% perturb = -2.55e7;
	T perturb = 0;
	//% same as  rf_pose_from_point_tangents_2_fn_t, but polynomial value might be
	//% perturbed by max of(perturb) so that roots appear for noisy data.This wasn't
	//% used in the end.
	//%

	//%function of t part :
	static T *fvalue = output[0];
	static T *A      = output[1];
	static T *B      = output[2];
	static T *C      = output[3];

	static T *E      = output[4];
	static T *F      = output[5];
	static T *G      = output[6];
	static T *H      = output[7];

	static T *J      = output[8];
	static T *K      = output[9];
	static T *L      = output[10];

	// `t` integer powers
	static T t_pow2 = t * t;
	static T t_pow3 = t_pow2 * t;
	static T t_pow4 = t_pow3 * t;
	static T t_pow5 = t_pow4 * t;
	static T t_pow6 = t_pow5 * t;
	static T t_pow7 = t_pow6 * t;
	static T t_pow8 = t_pow7 * t;

	// `(t^2 + 1)` integer powers
	static T t_pow2_plus1_pow2 = (t_pow2 + 1) * (t_pow2 + 1);
	static T t_pow2_plus1_pow3 = t_pow2_plus1_pow2 * (t_pow2 + 1);
	static T t_pow2_plus1_pow4 = t_pow2_plus1_pow3 * (t_pow2 + 1);

	// Denominators
	static T *A_den = t_pow2_plus1_pow2;
	static T *B_den = t_pow2_plus1_pow3;
	static T *C_den = t_pow2_plus1_pow4;
	static T *E_den = t_pow2_plus1_pow2;
	static T *F_den = t_pow2_plus1_pow3;
	static T *G_den = t_pow2_plus1_pow4;
	static T *H_den = t_pow2_plus1_pow4;
	static T *J_den = t_pow2_plus1_pow3;
	static T *K_den = t_pow2_plus1_pow3;
	static T *L_den = t_pow2_plus1_pow2;

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
	//       2nd power          3rd power               4th power
	static T A_pow2 = A * A,    A_pow3 = A_pow2 * A,    A_pow4 = A_pow3 * A;
	static T B_pow2 = B * B,    B_pow3 = B_pow2 * B,    B_pow4 = B_pow3 * B;
	static T C_pow2 = C * C,    C_pow3 = C_pow2 * C,    C_pow4 = C_pow3 * C;
	static T E_pow2 = E * E,    E_pow3 = E_pow2 * E,    E_pow4 = E_pow3 * E;
	static T F_pow2 = F * F,    F_pow3 = F_pow2 * F,    F_pow4 = F_pow3 * F;
	static T G_pow2 = G * G,    G_pow3 = G_pow2 * G,    G_pow4 = G_pow3 * G;
	static T H_pow2 = H * H,    H_pow3 = H_pow2 * H,    H_pow4 = H_pow3 * H;
	static T J_pow2 = J * J,    J_pow3 = J_pow2 * J,    J_pow4 = J_pow3 * J;
	static T K_pow2 = K * K,    K_pow3 = K_pow2 * K,    K_pow4 = K_pow3 * K;
	static T L_pow2 = L * L,    L_pow3 = L_pow2 * L,    L_pow4 = L_pow3 * L;

	static T fvalue_terms[52];

	fvalue_terms[0]  = E * E * B * B * H * H * J * J;
	fvalue_terms[1]  = G * G * C * C * L. ^ 4;
	fvalue_terms[2]  = G * G * A * A * K. ^ 4;
	fvalue_terms[3]  = E * E * A * A * H. ^ 4;
	fvalue_terms[4]  = E * E * C * C * J. ^ 4;
	fvalue_terms[5]  = -2 * E * A * H * H * G * C * L * L;
	fvalue_terms[6]  = 2 * E * E * A * H * H * C * J * J;
	fvalue_terms[7]  = -2 * E * E * C * J. ^ 3 * B * H;
	fvalue_terms[8]  = 2 * E * C * C * J * J * G * L * L;
	fvalue_terms[9]  = 2 * E * A * A * H * H * G * K * K;
	fvalue_terms[10] = -2 * E * E * A * H. ^ 3 * B * J;
	fvalue_terms[11] = -2 * E * A * H * H * G * B * K * L;
	fvalue_terms[12] = -2 * E * C * J * J * G * B * K * L;
	fvalue_terms[13] = -2 * E * C * J * J * G * A * K * K;
	fvalue_terms[14] = -2 * E * B * H * J * G * C * L * L;
	fvalue_terms[15] = -2 * E * B * H * J * G * A * K * K;
	fvalue_terms[16] = G * G * B * B * K * K * L * L;
	fvalue_terms[17] = -2 * G * G * B * K * L. ^ 3 * C;
	fvalue_terms[18] = -2 * G * G * B * K. ^ 3 * L * A;
	fvalue_terms[19] = 2 * G * G * C * L * L * A * K * K;
	fvalue_terms[20] = -2 * F * E * A * A * H. ^ 3 * K;
	fvalue_terms[21] = -2 * F * E * A * H * K * C * J * J;
	fvalue_terms[22] = 3 * F * E * A * H * H * K * B * J;
	fvalue_terms[23] = 3 * F * A * H * K * K * G * B * L;
	fvalue_terms[24] = -2 * F * A * H * K * G * C * L * L;
	fvalue_terms[25] = -2 * F * A * A * H * K. ^ 3 * G;
	fvalue_terms[26] = F * E * B * H. ^ 3 * L * A;
	fvalue_terms[27] = 3 * F * E * B * H * L * C * J * J;
	fvalue_terms[28] = -F * E * B * B * H * H * L * J;
	fvalue_terms[29] = -F * B * B * H * L * L * G * K;
	fvalue_terms[30] = F * B * H * L. ^ 3 * G * C;
	fvalue_terms[31] = F * E * B * K * J. ^ 3 * C;
	fvalue_terms[32] = -F * E * B * B * K * J * J * H;
	fvalue_terms[33] = -F * B * B * K * K * J * G * L;
	fvalue_terms[34] = 3 * F * B * K * J * G * C * L * L;
	fvalue_terms[35] = F * B * K. ^ 3 * J * G * A;
	fvalue_terms[36] = -2 * F * E * C * J * L * A * H * H;
	fvalue_terms[37] = -2 * F * E * C * C * J. ^ 3 * L;
	fvalue_terms[38] = -2 * F * C * C * J * L. ^ 3 * G;
	fvalue_terms[39] = -2 * F * C * J * L * G * A * K * K;
	fvalue_terms[40] = F * F * A * A * K * K * H * H;
	fvalue_terms[41] = F * F * A * K * K * C * J * J;
	fvalue_terms[42] = -F * F * A * K * K * B * H * J;
	fvalue_terms[43] = -F * F * B * K * L * A * H * H;
	fvalue_terms[44] = -F * F * B * K * L * C * J * J;
	fvalue_terms[45] = F * F * B * B * K * L * H * J;
	fvalue_terms[46] = F * F * C * L * L * A * H * H;
	fvalue_terms[47] = F * F * C * C * L * L * J * J;
	fvalue_terms[48] = -F * F * C * L * L * B * H * J;
	fvalue_terms[49] = G * E * B * B * H * H * L * L;
	fvalue_terms[50] = G * E * B * B * K * K * J * J;
	fvalue_terms[51] = 8 * G * E * A * H * K * C * J * L;

	for (int i = 0; i < 52; i++)
		fvalue += fvalue_terms[i];
}

}


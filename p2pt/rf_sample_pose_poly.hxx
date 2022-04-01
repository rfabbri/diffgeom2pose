#include <common.hxx>

//% samples/evaluates pose polynomial; 
template<typename T>
void rf_sample_pose_poly(
	const T t, 
	T A[2001], T B[2001], T C[2001], T E[2001], 
	T F[2001], T G[2001], T H[2001], T K[2001], 
	T J[2001], T K[2001], T L[2001], T fvalue[2001]
)
{
	// TODO: Check if `extern` is needed or another approach should be used
	//% function of t part:
	extern A0 A1 A2;
	extern B0 B1 B2 B3;
	extern C0 C1 C2 C3 C4;
	extern E0 E1 E2;
	extern F0 F1 F2 F3;
	extern G0 G1 G2 G3 G4;
	extern H0 H1 H2 H3 H4;
	extern J0 J1 J2 J3;
	extern K0 K1 K2 K3;
	extern L0 L1 L2;

	// Element-wise power t.^2, t.^3, t.^4, ...
	static T t_pow2[2001]; common::vec_el_wise_pow(t, 2, t_pow2);
	static T t_pow3[2001]; common::vec_el_wise_pow(t, 3, t_pow3);
	static T t_pow4[2001]; common::vec_el_wise_pow(t, 4, t_pow4);
	static T t_pow5[2001]; common::vec_el_wise_pow(t, 5, t_pow5);
	static T t_pow6[2001]; common::vec_el_wise_pow(t, 6, t_pow6);
	static T t_pow7[2001]; common::vec_el_wise_pow(t, 7, t_pow7);
	static T t_pow8[2001]; common::vec_el_wise_pow(t, 8, t_pow8);

	// Element-wise scalar addition (1 + t.^2)
	static T t_pow2_plus1[2001]; common::vec_add_scalar(t_pow2, 1, t_pow2_plus1);

	// Element-wise power (1 + t.^2)^2, (1 + t.^2)^3, (1 + t.^2)^4, ...
	static T t_pow2_plus1_pow2[2001]; common::vec_el_wise_pow(t_pow2_plus1, 2, t_pow2_plus1_pow2);
	static T t_pow2_plus1_pow3[2001]; common::vec_el_wise_pow(t_pow2_plus1, 3, t_pow2_plus1_pow3);
	static T t_pow2_plus1_pow4[2001]; common::vec_el_wise_pow(t_pow2_plus1, 4, t_pow2_plus1_pow4);

	// Denominators
	static T const *A_den = t_pow2_plus1_pow2;
	static T const *B_den = t_pow2_plus1_pow3;
	static T const *C_den = t_pow2_plus1_pow4;
	static T const *E_den = t_pow2_plus1_pow2;
	static T const *F_den = t_pow2_plus1_pow3;
	static T const *G_den = t_pow2_plus1_pow4;
	static T const *H_den = t_pow2_plus1_pow4;
	static T const *J_den = t_pow2_plus1_pow3;
	static T const *K_den = t_pow2_plus1_pow3;
	static T const *L_den = t_pow2_plus1_pow2;

	// Calculations
	// A
	for (int i = 0; i < 2001; i++) {
		A[i] = (A0 + A1 * t[i] + A2 * t_pow2[i] - A1 * t_pow3[i] + A0 * t_pow4[i]);
	}
	common::vec1vec2_el_wise_right_div(A, A_den, A);

	// B
	for (int i = 0; i < 2001; i++) {
		B[i] = (B0 + B1 * t[i] + B2 * t_pow2[i] + B3 * t_pow3[i] - B2 * t_pow4[i] + B1 * t_pow5[i] - B0 * t_pow6[i]);
	}
	common::vec1vec2_el_wise_right_div(B, B_den, B);

	// C
	for (int i = 0; i < 2001; i++) {
		C[i] = (C0 + C1 * t[i] + C2 * t_pow2[i] + C3 * t_pow3[i] + C4 * t_pow4[i] - C3 * t_pow5[i] + C2 * t_pow6[i] - C1 * t_pow7[i] + C0 * t_pow8[i]);
	}
	common::vec1vec2_el_wise_right_div(C, C_den, C);

	// E
	for (int i = 0; i < 2001; i++) {
		E[i] = (E0 + E1 * t[i] + E2 * t_pow2[i] - E1 * t_pow3[i] + E0 * t_pow4[i]);
	}
	common::vec1vec2_el_wise_right_div(E, E_den, E);

	// F
	for (int i = 0; i < 2001; i++) {
		F[i] = (F0 + F1 * t[i] + F2 * t_pow2[i] + F3 * t_pow3[i] - F2 * t_pow4[i] + F1 * t_pow5[i] - F0 * t_pow6[i]);
	}
	common::vec1vec2_el_wise_right_div(F, F_den, F);

	// G
	for (int i = 0; i < 2001; i++) {
		G[i] = (G0 + G1 * t[i] + G2 * t_pow2[i] + G3 * t_pow3[i] + G4 * t_pow4[i] - G3 * t_pow5[i] + G2 * t_pow6[i] - G1 * t_pow7[i] + G0 * t_pow8[i]);
	}
	common::vec1vec2_el_wise_right_div(G, G_den, G);

	// H
	for (int i = 0; i < 2001; i++) {
		H[i] = (H0 + H1 * t[i] + H2 * t_pow2[i] + H3 * t_pow3[i] + H4 * t_pow4[i] - H3 * t_pow5[i] + H2 * t_pow6[i] - H1 * t_pow7[i] + H0 * t_pow8[i]);
	}
	common::vec1vec2_el_wise_right_div(H, H_den, H);

	// J
	for (int i = 0; i < 2001; i++) {
		J[i] = (J0 + J1 * t[i] + J2 * t_pow2[i] + J3 * t_pow3[i] - J2 * t_pow4[i] + J1 * t_pow5[i] - J0 * t_pow6[i]);
	}
	common::vec1vec2_el_wise_right_div(J, J_den, J);

	// K
	for (int i = 0; i < 2001; i++) {
		K[i] = (K0 + K1 * t[i] + K2 * t_pow2[i] + K3 * t_pow3[i] - K2 * t_pow4[i] + K1 * t_pow5[i] - K0 * t_pow6[i]);
	}
	common::vec1vec2_el_wise_right_div(K, K_den, K);

	// L
	for (int i = 0; i < 2001; i++) {
		L[i] = (L0 + L1 * t[i] + L2 * t_pow2[i] - L1 * t_pow3[i] + L0 * t_pow4[i]);
	}
	common::vec1vec2_el_wise_right_div(L, L_den, L);

	static T A_pow2[2001]; common::vec_el_wise_pow(A, 2, A_pow2);
	static T B_pow2[2001]; common::vec_el_wise_pow(B, 2, B_pow2);
	static T C_pow2[2001]; common::vec_el_wise_pow(C, 2, C_pow2);
	static T E_pow2[2001]; common::vec_el_wise_pow(E, 2, E_pow2);
	static T F_pow2[2001]; common::vec_el_wise_pow(F, 2, F_pow2);
	static T G_pow2[2001]; common::vec_el_wise_pow(G, 2, G_pow2);
	static T H_pow2[2001]; common::vec_el_wise_pow(H, 2, H_pow2);
	static T J_pow2[2001]; common::vec_el_wise_pow(J, 2, J_pow2);
	static T K_pow2[2001]; common::vec_el_wise_pow(K, 2, K_pow2);
	static T L_pow2[2001]; common::vec_el_wise_pow(L, 2, L_pow2);

	static T H_pow3[2001]; common::vec_el_wise_pow(H, 3, H_pow3);
	static T J_pow3[2001]; common::vec_el_wise_pow(J, 3, J_pow3);
	static T K_pow3[2001]; common::vec_el_wise_pow(K, 3, K_pow3);
	static T L_pow3[2001]; common::vec_el_wise_pow(L, 3, L_pow3);

	static T H_pow4[2001]; common::vec_el_wise_pow(H, 4, H_pow4);
	static T J_pow4[2001]; common::vec_el_wise_pow(J, 4, J_pow4);
	static T K_pow4[2001]; common::vec_el_wise_pow(K, 4, K_pow4);
	static T L_pow4[2001]; common::vec_el_wise_pow(L, 4, L_pow4);

	// `fvalue` is composed of the sum of these terms
	static T fvalue_terms[52][2001];

	// fvalue_terms[0] 
	common::vec1vec2_el_wise_mult(E_pow2,
	common::vec1vec2_el_wise_mult(B_pow2, 
	common::vec1vec2_el_wise_mult(H_pow2, J_pow2, fvalue_terms[0]), fvalue_terms[0]), fvalue_terms[0]);

	// fvalue_terms[1]
	common::vec1vec2_el_wise_mult(G_pow2,
	common::vec1vec2_el_wise_mult(C_pow2, L_pow4, fvalue_terms[1]), fvalue_terms[1]);

	// fvalue_terms[2]
	common::vec1vec2_el_wise_mult(G_pow2,
	common::vec1vec2_el_wise_mult(A_pow2, K_pow4, fvalue_terms[2]), fvalue_terms[2]);

	// fvalue_terms[3]
	common::vec1vec2_el_wise_mult(E_pow2,
	common::vec1vec2_el_wise_mult(A_pow2, H_pow4, fvalue_terms[3]), fvalue_terms[3]);

	// fvalue_terms[4]
	common::vec1vec2_el_wise_mult(E_pow2,
	common::vec1vec2_el_wise_mult(C_pow2, J_pow4, fvalue_terms[4]), fvalue_terms[4]);

	// fvalue_terms[5]
	// fvalue_terms[6]
	// fvalue_terms[7]
	// fvalue_terms[8]
	// fvalue_terms[9]
	// fvalue_terms[10]
	// fvalue_terms[11]
	// fvalue_terms[12]
	// fvalue_terms[13]
	// fvalue_terms[14]
	// fvalue_terms[15]
	// fvalue_terms[16]
	// fvalue_terms[17]
	// fvalue_terms[18]
	// fvalue_terms[19]
	// fvalue_terms[20]
	// fvalue_terms[21]
	// fvalue_terms[22]
	// fvalue_terms[23]
	// fvalue_terms[24]
	// fvalue_terms[25]
	// fvalue_terms[26]
	// fvalue_terms[27]
	// fvalue_terms[28]
	// fvalue_terms[29]
	// fvalue_terms[30]
	// fvalue_terms[31]
	// fvalue_terms[32]
	// fvalue_terms[33]
	// fvalue_terms[34]
	// fvalue_terms[35]
	// fvalue_terms[36]
	// fvalue_terms[37]
	// fvalue_terms[38]
	// fvalue_terms[39]
	// fvalue_terms[40]
	// fvalue_terms[41]
	// fvalue_terms[42]
	// fvalue_terms[43]
	// fvalue_terms[44]
	// fvalue_terms[45]
	// fvalue_terms[46]
	// fvalue_terms[47]
	// fvalue_terms[48]
	// fvalue_terms[49]
	// fvalue_terms[50]
	// fvalue_terms[51]

	for (int i = 0; i < 52; i++) {
		for (int j = 0; j < 2001; j++) {
			fvalue[j] += fvalue_terms[i][j];
		}
	}

	// TODO: rewrite `fvalue`
	fvalue = 
		/*  0 */ E .* E .* B .* B .* H .* H .* J .* J 
		/*  1 */ + G .* G .* C .* C .* L.^4 
		/*  2 */ + G .* G .* A .* A .* K.^4 
		/*  3 */ + E .* E .* A .* A .* H.^4 
		/*  4 */ + E .* E .* C .* C .* J.^4
		/*  5 */ - 2 .* E .* A .* H .* H .* G .* C .* L .* L 
		/*  6 */ + 2 .* E .* E .* A .* H .* H .* C .* J .* J 
		/*  7 */ - 2 .* E .* E .* C .* J.^3 .* B .* H 
		/*  8 */ + 2 .* E .* C .* C .* J .* J .* G .* L .* L 
		/*  9 */ + 2 .* E .* A .* A .* H .* H .* G .* K .* K 
		/* 10 */ - 2 .* E .* E .* A .* H.^3 .* B .* J 
		/* 11 */ - 2 .* E .* A .* H .* H .* G .* B .* K .* L 
		/* 12 */ - 2 .* E .* C .* J .* J .* G .* B .* K .* L 
		/* 13 */ - 2 .* E .* C .* J .* J .* G .* A .* K .* K 
		/* 14 */ - 2 .* E .* B .* H .* J .* G .* C .* L .* L 
		/* 15 */ - 2 .* E .* B .* H .* J .* G .* A .* K .* K 
		/* 16 */ + G .* G .* B .* B .* K .* K .* L .* L 
		/* 17 */ - 2 .* G .* G .* B .* K .* L.^3 .* C 
		/* 18 */ - 2 .* G .* G .* B .* K.^3 .* L .* A 
		/* 19 */ + 2 .* G .* G .* C .* L .* L .* A .* K .* K 
		/* 20 */ - 2 .* F .* E .* A .* A .* H.^3 .* K 
		/* 21 */ - 2 .* F .* E .* A .* H .* K .* C .* J .* J 
		/* 22 */ + 3 .* F .* E .* A .* H .* H .* K .* B .* J 
		/* 23 */ + 3 .* F .* A .* H .* K .* K .* G .* B .* L 
		/* 24 */ - 2 .* F .* A .* H .* K .* G .* C .* L .* L 
		/* 25 */ - 2 .* F .* A .* A .* H .* K.^3 .* G 
		/* 26 */ + F .* E .* B .* H.^3 .* L .* A 
		/* 27 */ + 3 .* F .* E .* B .* H .* L .* C .* J .* J 
		/* 28 */ - F .* E .* B .* B .* H .* H .* L .* J 
		/* 29 */ - F .* B .* B .* H .* L .* L .* G .* K 
		/* 30 */ + F .* B .* H .* L.^3 .* G .* C 
		/* 31 */ + F .* E .* B .* K .* J.^3 .* C 
		/* 32 */ - F .* E .* B .* B .* K .* J .* J .* H 
		/* 33 */ - F .* B .* B .* K .* K .* J .* G .* L 
		/* 34 */ + 3 .* F .* B .* K .* J .* G .* C .* L .* L 
		/* 35 */ + F .* B .* K.^3 .* J .* G .* A 
		/* 36 */ - 2 .* F .* E .* C .* J .* L .* A .* H .* H 
		/* 37 */ - 2 .* F .* E .* C .* C .* J.^3 .* L 
		/* 38 */ - 2 .* F .* C .* C .* J .* L.^3 .* G 
		/* 39 */ - 2 .* F .* C .* J .* L .* G .* A .* K .* K 
		/* 40 */ + F .* F .* A .* A .* K .* K .* H .* H 
		/* 41 */ + F .* F .* A .* K .* K .* C .* J .* J 
		/* 42 */ - F .* F .* A .* K .* K .* B .* H .* J 
		/* 43 */ - F .* F .* B .* K .* L .* A .* H .* H 
		/* 44 */ - F .* F .* B .* K .* L .* C .* J .* J 
		/* 45 */ + F .* F .* B .* B .* K .* L .* H .* J 
		/* 46 */ + F .* F .* C .* L .* L .* A .* H .* H 
		/* 47 */ + F .* F .* C .* C .* L .* L .* J .* J 
		/* 48 */ - F .* F .* C .* L .* L .* B .* H .* J 
		/* 49 */ + G .* E .* B .* B .* H .* H .* L .* L 
		/* 50 */ + G .* E .* B .* B .* K .* K .* J .* J 
		/* 51 */ + 8 .* G .* E .* A .* H .* K .* C .* J .* L;

}
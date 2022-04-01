#include <common.hxx>

//% samples/evaluates pose polynomial; 
template<typename T>
void rf_sample_pose_poly(
	const T t, 
	const T A[2001], const T B[2001], const T C[2001], const T D[2001],
	const T E[2001], const T F[2001], const T G[2001], const T H[2001], 
	const T K[2001], const T L[2001], const T M[2001], const T fvalue[2001]
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


	const double t_pow2[2001];
	const double t_pow3[2001];
	const double t_pow4[2001];
	const double t_pow5[2001];
	const double t_pow6[2001];
	const double t_pow7[2001];
	const double t_pow8[2001];
	common::vec_el_wise_pow(t, 2, t_pow2);
	common::vec_el_wise_pow(t, 3, t_pow3);
	common::vec_el_wise_pow(t, 4, t_pow4);
	common::vec_el_wise_pow(t, 5, t_pow5);
	common::vec_el_wise_pow(t, 6, t_pow6);
	common::vec_el_wise_pow(t, 7, t_pow7);
	common::vec_el_wise_pow(t, 8, t_pow8);

	const double t_pow2_plus1[2001]; 
	common::vec_add_scalar(t_pow2, 1, t_pow2_plus1);

	const double t_pow2_plus1_pow2[2001];
	const double t_pow2_plus1_pow3[2001];
	const double t_pow2_plus1_pow4[2001];
	common::vec_el_wise_pow(t_pow2_plus1, 2, t_pow2_plus1_pow2);
	common::vec_el_wise_pow(t_pow2_plus1, 3, t_pow2_plus1_pow3);
	common::vec_el_wise_pow(t_pow2_plus1, 4, t_pow2_plus1_pow4);

	// A
	for (int i = 0; i < 2001; i++) {
		A[i] = (A0 + A1 * t[i] + A2 * t_pow2[i] - A1 * t_pow3[i] + A0 * t_pow4[i]);
	}
	common::vec1vec2_el_wise_right_div(A, t_pow2_plus1_pow2, A);

	// B
	for (int i = 0; i < 2001; i++) {
		B[i] = (B0 + B1 * t[i] + B2 * t_pow2[i] + B3 * t_pow3[i] - B2 * t_pow4[i] + B1 * t_pow5[i] - B0 * t_pow6[i]);
	}
	common::vec1vec2_el_wise_right_div(B, t_pow2_plus1_pow3, B);

	// TODO: Modify MATLAB operators `.^` and `.*`
	A = (A0 + A1*t + A2*t.^2 - A1*t.^3 + A0*t.^4) ./ (1+t.^2).^2;
	B = (B0 + B1 *t + B2 *t.^2 + B3 *t.^3 - B2 *t.^4 + B1 *t.^5 - B0 *t.^6) ./ (1+t.^2).^3;
	C = (C0 + C1 *t + C2 *t.^2 + C3 *t.^3 + C4 *t.^4 - C3 *t.^5 + C2 *t.^6 - C1 *t.^7 + C0 *t.^8)./(1+t.^2).^4;
	E = (E0 + E1 *t + E2 *t.^2 - E1 *t.^3 + E0 *t.^4)./(1+t.^2).^2;
	F = (F0 + F1 *t + F2 *t.^2 + F3 *t.^3 - F2 *t.^4 + F1 *t.^5 - F0 *t.^6)./(1+t.^2).^3;
	G = (G0 + G1 *t + G2 *t.^2 + G3 *t.^3 + G4 *t.^4 - G3 *t.^5 + G2 *t.^6 - G1 *t.^7 + G0 *t.^8) ./ (1+t.^2).^4;
	H = (H0 + H1 *t + H2 *t.^2 + H3 *t.^3 + H4 *t.^4 - H3 *t.^5 + H2 *t.^6 - H1 *t.^7 + H0 *t.^8)./(1+t.^2).^4;
	J = (J0 + J1 *t + J2 *t.^2 + J3 *t.^3 - J2 *t.^4 + J1 *t.^5 - J0 *t.^6) ./ (1+t.^2).^3 ;
	K = (K0 + K1 *t + K2 *t.^2 + K3 *t.^3 - K2 *t.^4 + K1 *t.^5 - K0 *t.^6) ./ (1+t.^2).^3;
	L = (L0 + L1 *t + L2 *t.^2 - L1 *t.^3 + L0 *t.^4) ./ (1+t.^2).^2;

	fvalue = E .* E .* B .* B .* H .* H .* J .* J + G .* G .* C .* C .* L.^4 + G .* G .* A ...
	.* A .* K.^4 + E .* E .* A .* A .* H.^4 + E .* E .* C .* C .* J.^4 ...
	- 2 .* E .* A .* H .* H .* G .* C .* L .* L + 2 .* E .* E .* A .* H .* H .* C .* ...
	J .* J - 2 .* E .* E .* C .* J.^3 .* B .* H + 2 .* E .* C .* C .* J .* J .* ...
	G .* L .* L + 2 .* E .* A .* A .* H .* H .* G .* K .* K - 2 .* E .* E .* ...
	A .* H.^3 .* B .* J - 2 .* E .* A .* H .* H .* G .* B .* K .* L - 2 .* E .* C .* J .* J .* ...
	G .* B .* K .* L - 2 .* E .* C .* J .* J .* G .* A .* K .* K - 2 .* E .* B .* H .* J .* ...
	G .* C .* L .* L - 2 .* E .* B .* H .* J .* G .* A .* K .* K + G .* G .* B .* B .* K .* K .* ...
	L .* L - 2 .* G .* G .* B .* K .* L.^3 .* C - 2 .* G .* G .* B .* K.^3 .* L .* A + 2 .* ...
	G .* G .* C .* L .* L .* A .* K .* K - 2 .* F .* E .* A .* A .* ...
	H.^3 .* K - 2 .* F .* E .* A .* H .* K .* C .* J .* J + 3 .* F .* E .* A .* ...
	H .* H .* K .* B .* J + 3 .* F .* A .* H .* K .* K .* G .* B .* L - 2 .* F .* A .* H .* ...
	K .* G .* C .* L .* L - 2 .* F .* A .* A .* H .* K.^3 .* G + F .* E .* B .* ...
	H.^3 .* L .* A + 3 .* F .* E .* B .* H .* L .* C .* J .* J - F .* E .* B .* B .* ...
	H .* H .* L .* J - F .* B .* B .* H .* L .* L .* G .* K + F .* B .* H .* L.^3 .* G .* ...
	C + F .* E .* B .* K .* J.^3 .* C - F .* E .* B .* B .* K .* J .* J .* H - F .* B .* ...
	B .* K .* K .* J .* G .* L + 3 .* F .* B .* K .* J .* G .* C .* L .* L + F...
	.* B .* K.^3 .* J .* G .* A - 2 .* F .* E .* C .* J .* L .* A .* H .* H - 2 .* F .* E .* C .* ...
	C .* J.^3 .* L - 2 .* F .* C .* C .* J .* L.^3 .* G - 2 .* F .* ...
	C .* J .* L .* G .* A .* K .* K + F .* F .* A .* A .* K .* K .* H .* H + F .* F .* A .* K .* K .* ...
	C .* J .* J - F .* F .* A .* K .* K .* B .* H .* J - F .* F .* B .* K .* L .* A .* H .* H - F .* ...
	F .* B .* K .* L .* C .* J .* J + F .* F .* B .* B .* K .* L .* H .* J + F .* F .* C .* L .* L .* ...
	A .* H .* H + F .* F .* C .* C .* L .* L .* J .* J - F .* F .* C .* L .* L .* B .* H .* J + G .* ...
	E .* B .* B .* H .* H .* L .* L + G .* E .* B .* B .* K .* K .* J .* J + 8 .* G .* E .* A .* ...
	H .* K .* C .* J .* L;

}
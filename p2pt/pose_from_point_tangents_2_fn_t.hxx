#ifndef pose_from_point_tangents_2_fn_t_hxx_
#define pose_from_point_tangents_2_fn_t_hxx_

namespace P2Pt {

template<typename T>
inline T
pose_poly<T>::
pose_from_point_tangents_2_fn_t(const T t, T (*output)[10] /* = nullptr */)
{
	T buf[10];

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

	const T t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t, t6 = t5*t, t7 = t6*t, t8 = t7*t;

	// `(t^2 + 1)` integer powers
  // TODO(check) when removing static, the tests accuse loss of precision
	const T 
    t2_plus1_pow2 = (t2 + 1.) * (t2 + 1.),
    t2_plus1_pow3 = t2_plus1_pow2 * (t2 + 1.),
    t2_plus1_pow4 = t2_plus1_pow3 * (t2 + 1.);

	A = (A0+A1*t+A2*t2-A1*t3+A0*t4)/t2_plus1_pow2;
	B = (B0+B1*t+B2*t2+B3*t3-B2*t4+B1*t5-B0*t6)/t2_plus1_pow3;
	C = (C0+C1*t+C2*t2+C3*t3+C4*t4-C3*t5+C2*t6-C1*t7+C0*t8)/t2_plus1_pow4;
	E = (E0+E1*t+E2*t2-E1*t3+E0*t4)/t2_plus1_pow2;
	F = (F0+F1*t+F2*t2+F3*t3-F2*t4+F1*t5-F0*t6)/t2_plus1_pow3;
	G = (G0+G1*t+G2*t2+G3*t3+G4*t4-G3*t5+G2*t6-G1*t7+G0*t8)/t2_plus1_pow4;
	H = (H0+H1*t+H2*t2+H3*t3+H4*t4-H3*t5+H2*t6-H1*t7+H0*t8)/t2_plus1_pow4;
	J = (J0+J1*t+J2*t2+J3*t3-J2*t4+J1*t5-J0*t6)/t2_plus1_pow3;
	K = (K0+K1*t+K2*t2+K3*t3-K2*t4+K1*t5-K0*t6)/t2_plus1_pow3;
	L = (L0+L1*t+L2*t2-L1*t3+L0*t4)/t2_plus1_pow2;

	const T
	A2 = A * A,
	B2 = B * B,
	C2 = C * C,
	E2 = E * E,
	F2 = F * F,
	G2 = G * G,
	H2 = H * H,    H3 = H2 * H,    H4 = H3 * H,
	J2 = J * J,    J3 = J2 * J,   
	K2 = K * K,    K3 = K2 * K,
  L2 = L*L,      L3 = L2 * L;

  return // fvalue_terms[X]
	// TODO: Analyze the need of extra/lower precison for the calculation of these terms.
	/*0*/E2*B2*H2*J2
	/*1*/+G2*C2*L*L*L*L
	/*2*/+G2*A2*K3*K
	/*3*/+E2*A2*H4
	/*4*/+E2*C2*J3*J
	/*5*/+-2.*E*A*H2*G*C*L2
	/*6*/+2.*E2*A*H2*C*J2
	/*7*/+-2.*E2*C*J3*B*H
	/*8*/+2.*E*C2*J2*G*L2
	/*9*/+2.*E*A2*H2*G*K2
	/*10*/+-2.*E2*A*H3*B*J
	/*11*/+-2.*E*A*H2*G*B*K*L
	/*12*/+-2.*E*C*J2*G*B*K*L
	/*13*/+-2.*E*C*J2*G*A*K2
	/*14*/+-2.*E*B*H*J*G*C*L2
	/*15*/+-2.*E*B*H*J*G*A*K2
	/*16*/+G2*B2*K2*L2
	/*17*/+-2.*G2*B*K*L3*C
	/*18*/+-2.*G2*B*K3*L*A
	/*19*/+2.*G2*C*L2*A*K2
	/*20*/+-2.*F*E*A2*H3*K
	/*21*/+-2.*F*E*A*H*K*C*J2
	/*22*/+3.*F*E*A*H2*K*B*J
	/*23*/+3.*F*A*H*K2*G*B*L
	/*24*/+-2.*F*A*H*K*G*C*L2
	/*25*/+-2.*F*A2*H*K3*G
	/*26*/+F*E*B*H3*L*A
	/*27*/+3.*F*E*B*H*L*C*J2
	/*28*/+-1.*F*E*B2*H2*L*J
	/*29*/+-1.*F*B2*H*L2*G*K
	/*30*/+F*B*H*L3*G*C
	/*31*/+F*E*B*K*J3*C
	/*32*/+-1.*F*E*B2*K*J2*H
	/*33*/+-1.*F*B2*K2*J*G*L
	/*34*/+3.*F*B*K*J*G*C*L2
	/*35*/+F*B*K3*J*G*A
	/*36*/+-2.*F*E*C*J*L*A*H2
	/*37*/+-2.*F*E*C2*J3*L
	/*38*/+-2.*F*C2*J*L3*G
	/*39*/+-2.*F*C*J*L*G*A*K2
	/*40*/+F2*A2*K2*H2
	/*41*/+F2*A*K2*C*J2
	/*42*/+-1.*F2*A*K2*B*H*J
	/*43*/+-1.*F2*B*K*L*A*H2
	/*44*/+-1.*F2*B*K*L*C*J2
	/*45*/+F2*B2*K*L*H*J
	/*46*/+F2*C*L2*A*H2
	/*47*/+F2*C2*L2*J2
	/*48*/+-1.*F2*C*L2*B*H*J
	/*49*/+G*E*B2*H2*L2
	/*50*/+G*E*B2*K2*J2
	/*51*/+8.*G*E*A*H*K*C*J*L;
}

}

#endif // !pose_from_point_tangents_2_fn_t_hxx_


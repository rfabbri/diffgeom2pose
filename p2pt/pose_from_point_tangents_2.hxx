#ifndef pose_from_point_tangents_2_hxx_
#define pose_from_point_tangents_2_hxx_

namespace P2Pt {

using namespace common;

template<typename T>
void
pose_poly<T>::
pose_from_point_tangents_2(
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3]
)
{
	static constexpr T PI = 3.141592653589793;

	T g11, g12, g21, g22, h11, h12, h21, h22;
	static T V[3], buf[3];
	static T a1, a2, a3, a4, a5, a6;
	static T t1, t2, t4, t5, t6, t7, t8;
	static T t11, t14, /* t15, */ t21, t24, t25;
	static T den1, den2;

	g11 = gama1[0]; g12 = gama1[1];
	g21 = gama2[0]; g22 = gama2[1];

	h11 = tgt1[0];  h12 = tgt1[1];
	h21 = tgt2[0];  h22 = tgt2[1];

	vec1vec2_3el_sub(Gama1, Gama2, V);
	vec_3el_wise_mult2(V, V, buf);       a1 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt1, Tgt1, buf); a2 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt2, Tgt2, buf); a3 = vec_3el_sum(buf);
	vec_3el_wise_mult2(V, Tgt1, buf);    a4 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt1, Tgt2, buf); a5 = vec_3el_sum(buf);
	vec_3el_wise_mult2(V, Tgt2, buf);    a6 = vec_3el_sum(buf);

	t4 = g11 * g11;
	t5 = g12 * g12;
	t6 = g21 * g21;
	t7 = g22 * g22;
	t11 = 2 * (1 + g11 * g21 + g12 * g22) / (t4 + t5 - t6 - t7);

	theta = 0.5*atan(t11);
	if (theta < 0) theta += PI / 2;
	sth = sin(theta);
	cth = cos(theta);

	//% 497-798
	//%theta = .7865071740;
	//
	//% 240-1100
	//%theta = .7887953040;
	//
	//% 101-406:
	//%theta = .7852237735
	//
	//% the above theta has similar sin(2theta), cos(2theta) as the maple spreadsheet.

	t1 = sin(2*theta);
	t2 = cos(2*theta);
	t5 = g22 * g22;
	t6 = t2 * t2;
	t8 = g21 * g21;
	t14 = g12 * g12;
	//double t15 = t1 * t1; // double-checked: not used in matlab
	t21 = g11 * g11;

	den1 = 2*t1*(g11*g21 +g12*g22 + 1) + t2*(t21 + t14 - t8 - t5);
	den2 = t21 + t14 + t8 + t5 + 2;

	t25 = -2*a1 / (den1 - den2);
	beta = sqrt(t25);

	t24 = 2*a1 / (den1 + den2);
	alpha = sqrt(t24);

	//% Coefficient code adapted from Maple ::: can be further cleaned up but works

	A0 = a4 * a4 * g12 * g12
	+ a4 * a4 * g11 * g11
	+ a4 * a4
	+ 2.0 * a2 * intpow(g11, 3) * g21 * beta * beta * sth * cth
	+ 2.0 * a2 * g21 * g11 * g12 * g12 * beta * beta * sth * cth
	- 2.0 * a2 * g11 * g11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- a2 * intpow(g12, 4) * beta * beta * intpow(sth, 2)
	- a2 * g21 * g21 * g11 * g11 * beta * beta * intpow(cth, 2)
	+ 2.0 * a2 * g12 * g12 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * g22 * g12 * beta * beta * sth * cth
	- a2 * beta * beta * intpow(cth, 2)
	+ 2.0 * a2 * intpow(g12, 3) * g22 * beta * beta * sth * cth
	- a2 * intpow(g11, 4) * beta * beta * intpow(sth, 2)
	- 2.0 * a2 * g11 * g11 * beta * beta * intpow(sth, 2)
	- 2.0 * a2 * g12 * g12 * beta * beta * intpow(sth, 2)
	+ 2.0 * a2 * beta * beta * sth * cth
	- 2.0 * a2 * g21 * g11 * g22 * g12 * beta * beta * intpow(cth, 2)
	- a2 * beta * beta * intpow(sth, 2)
	+ 2.0 * a2 * g21 * g11 * beta * beta * sth * cth
	- a2 * g22 * g22 * g12 * g12 * beta * beta * intpow(cth, 2)
	- 2.0 * a2 * g22 * g12 * beta * beta * intpow(cth, 2)
	- 2.0 * a2 * g21 * g11 * beta * beta * intpow(cth, 2)
	+ 2.0 * a2 * g22 * g12 * beta * beta * sth * cth;

	A1 = 0.4e1 * a2 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * g21 * g11 * g11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.8e1 * a2 * g21 * g11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * g12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g22 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * g22 * g12 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g12 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * g11 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * g11 * g12 * g12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g21 * g11 * beta * intpow(sth, 2) * alpha
	- 0.8e1 * a2 * g11 * g11 * g12 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * intpow(g11, 4) * alpha * sth * beta * cth
	- 0.8e1 * a2 * g11 * g11 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * g11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g12, 3) * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * intpow(g12, 3) * g22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * intpow(g12, 4) * alpha * sth * beta * cth
	- 0.8e1 * a2 * g12 * g12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g11, 3) * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * intpow(g11, 3) * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * g11 * g22 * g12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * g22 * g12 * beta * intpow(sth, 2) * alpha;

	A2 = (2 * a4 * a4 * g12 * g12)
	+ (2 * a4 * a4 * g11 * g11)
	+ (2 * a4 * a4)
	+ 0.2e1 * a2 * intpow(g12, 4) * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * intpow(g11, 4) * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * (g11 * g11) * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * (g12 * g12) * beta * beta * intpow(sth, 2)
	- 0.4e1 * a2 * beta * beta * sth * cth
	+ 0.2e1 * a2 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * g21 * g11 * (g12 * g12) * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g21 * (g11 * g11) * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * (g12 * g12) * beta * beta * sth * cth
	+ 0.4e1 * a2 * g21 * g11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g22 * g22 * (g12 * g12) * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * g22 * g12 * beta * beta * sth * cth
	- 0.4e1 * a2 * (g11 * g11) * beta * beta * sth * cth
	- 0.4e1 * a2 * g21 * g11 * beta * beta * sth * cth
	+ 0.4e1 * a2 * (g11 * g11) * (g12 * g12) * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g22 * g12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * intpow(g11, 3) * g21 * beta * beta * sth * cth
	- 0.4e1 * a2 * (g11 * g11) * g22 * g12 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g21 * g11 * g22 * g12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * intpow(g12, 3) * g22 * beta * beta * sth * cth
	- 0.4e1 * a2 * intpow(g11, 4) * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a2 * (g11 * g11) * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a2 * intpow(g12, 4) * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a2 * (g12 * g12) * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a2 * alpha * alpha * cth * sth
	- 0.4e1 * a2 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a2 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a2 * g22 * g12 * alpha * alpha * cth * sth
	- 0.4e1 * a2 * g21 * g21 * (g11 * g11) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a2 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a2 * g22 * g22 * (g12 * g12) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a2 * g21 * g11 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * (g11 * g11) * (g12 * g12) * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a2 * (g11 * g11) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * g22 * g12 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a2 * intpow(g12, 3) * g22 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g22 * g12 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a2 * intpow(g11, 3) * g21 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * (g11 * g11) * g22 * g12 * alpha * alpha * cth * sth;

	B0 = -0.2e1 * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * intpow(g12, 3) * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * intpow(g11, 3) * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * intpow(cth, 2)
	+ a2 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a4 * a4 * h12 * g12);

	B1 = -0.2e1 * beta * sth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 0.4e1 * alpha * cth * (a2 * g21 * g11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * intpow(g12, 3) * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * intpow(g11, 3) * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * intpow(cth, 2)
	+ a2 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a4 * a4 * h12 * g12);

	B2 = -0.2e1 * beta * sth * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 0.2e1 * a4 * a4 * h11 * g11
	- 0.2e1 * a4 * a4 * h12 * g12
	- 0.2e1 * a2 * g22 * g22 * h12 * g12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * g11 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * h11 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * intpow(g11, 3) * h11 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * intpow(g12, 3) * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g21 * h11 * g22 * g12 * beta * beta * intpow(cth, 2))
	- 0.4e1 * alpha * cth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	+ 0.2e1 * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * intpow(g12, 3) * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * intpow(cth, 2)
	+ a2 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a2 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * intpow(g11, 3) * h11 * beta * beta * intpow(sth, 2)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * intpow(cth, 2)
	+ a2 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- a4 * a4 * h12 * g12);

	B3 = -0.2e1 * beta * sth * (-0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g12 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g21 * h11 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g21 * h11 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 0.4e1 * alpha * cth * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 0.2e1 * a4 * a4 * h11 * g11
	- 0.2e1 * a4 * a4 * h12 * g12
	- 0.2e1 * a2 * g22 * g22 * h12 * g12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * g12 * g12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * g11 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a2 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * h11 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * intpow(g11, 3) * h11 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * intpow(g12, 3) * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g21 * h11 * g22 * g12 * beta * beta * intpow(cth, 2))
	+ 0.2e1 * beta * sth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g12, 3) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * intpow(g11, 3) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth);

	C0 = -beta * beta * intpow(sth, 2) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2));

	C1 = -beta * beta * intpow(sth, 2) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- 0.4e1 * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2));

	C2 = -beta * beta * intpow(sth, 2) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2))
	- 0.4e1 * beta * sth * alpha * cth * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- (-0.2e1 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * alpha * alpha * intpow(cth, 2)) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2));

	C3 = -beta * beta * intpow(sth, 2) * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- 0.4e1 * beta * sth * alpha * cth * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2))
	- (-0.2e1 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * alpha * alpha * intpow(cth, 2)) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	+ 0.4e1 * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2));

	C4 = -0.2e1 * beta * beta * intpow(sth, 2) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2))
	- 0.4e1 * beta * sth * alpha * cth * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- (-0.2e1 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * alpha * alpha * intpow(cth, 2)) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * intpow(cth, 2))
	+ 0.4e1 * beta * sth * alpha * cth * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth);

	E0 = 0.2e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * intpow(g22, 3) * beta * beta * cth * sth
	- a3 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g12 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * g21 * beta * beta * intpow(sth, 2)
	- a3 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g11 * g21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * g22 * beta * beta * cth * sth
	- a3 * intpow(g21, 4) * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * g21 * beta * beta * cth * sth
	- a3 * g12 * g12 * g22 * g22 * beta * beta * intpow(sth, 2)
	- a3 * g11 * g11 * g21 * g21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ a6 * a6 * g21 * g21
	+ a6 * a6 * g22 * g22
	+ 0.2e1 * a3 * g11 * intpow(g21, 3) * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ a6 * a6
	- 0.2e1 * a3 * g21 * g21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a3 * intpow(g22, 4) * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * beta * beta * cth * sth;

	E1 = -0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.8e1 * a3 * g11 * g21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g12 * intpow(g22, 3) * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * intpow(g22, 3) * beta * intpow(cth, 2) * alpha
	+ 0.8e1 * a3 * g21 * g21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * intpow(g21, 4) * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * g21 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * g21 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g12 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * g21 * g22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * intpow(g21, 3) * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * intpow(g21, 3) * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.8e1 * a3 * g22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g21 * g21 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * intpow(g22, 4) * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * g21 * g22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * beta * intpow(cth, 2) * alpha;

	E2 = -0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	- 0.4e1 * a3 * g12 * g22 * beta * beta * cth * sth
	- 0.4e1 * a3 * g12 * intpow(g22, 3) * beta * beta * cth * sth
	+ 0.2e1 * a3 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a3 * g11 * g21 * beta * beta * cth * sth
	- 0.4e1 * a3 * g22 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * intpow(g21, 4) * beta * beta * intpow(cth, 2)
	- 0.4e1 * a3 * g21 * g21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * g12 * g22 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * g11 * g21 * g21 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a6 * a6 * g21 * g21
	+ 0.2e1 * a6 * a6 * g22 * g22
	- 0.4e1 * a3 * g11 * intpow(g21, 3) * beta * beta * cth * sth
	- 0.4e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a6 * a6
	- 0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a3 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a3 * g11 * g21 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a3 * g12 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g12 * intpow(g22, 3) * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a3 * g22 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g21 * g21 * g12 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g11 * intpow(g21, 3) * alpha * alpha * sth * cth
	- 0.4e1 * a3 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a3 * g21 * g21 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g11 * g21 * alpha * alpha * sth * cth
	- 0.4e1 * a3 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a3 * g11 * g21 * g22 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * intpow(g22, 4) * beta * beta * intpow(cth, 2)
	- 0.4e1 * a3 * beta * beta * cth * sth
	- 0.4e1 * a3 * intpow(g21, 4) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a3 * g21 * g21 * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a3 * alpha * alpha * sth * cth
	- 0.4e1 * a3 * intpow(g22, 4) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a3 * g22 * g22 * alpha * alpha * intpow(sth, 2);

	F0 = -0.2e1 * beta * cth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * intpow(sth, 2)
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * intpow(sth, 2)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * intpow(sth, 2)
	+ a3 * intpow(g21, 3) * h21 * beta * beta * intpow(cth, 2)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * intpow(sth, 2));

	F1 = -0.2e1 * beta * cth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * intpow(cth, 2) * alpha)
	+ 0.4e1 * alpha * sth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * intpow(sth, 2)
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * intpow(sth, 2)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * intpow(sth, 2)
	+ a3 * intpow(g21, 3) * h21 * beta * beta * intpow(cth, 2)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * intpow(sth, 2));

	F2 = -0.2e1 * beta * cth * (-(2 * a6 * a6 * h22 * g22)
	- (2 * a6 * a6 * h21 * g21)
	+ 0.2e1 * a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	- 0.2e1 * a3 * intpow(g22, 3) * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g21 * h21 * (g22 * g22) * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * h21 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * (g21 * g21) * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * h21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * intpow(g21, 3) * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * (g21 * g21) * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * (g21 * g21) * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * (g22 * g22) * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g11 * h21 * (g21 * g21) * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * (g22 * g22) * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g12 * h22 * (g22 * g22) * alpha * alpha * sth * cth)
	+ 0.4e1 * alpha * sth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * (g21 * g21) * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * (g22 * g22) * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * (g22 * g22) * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g21 * h21 * (g22 * g22) * alpha * sth * beta * cth
	+ 0.2e1 * a3 * (g21 * g21) * g12 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * (g21 * g21) * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * (g21 * g21) * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * (g22 * g22) * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * (g21 * g21) * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g11 * h21 * (g22 * g22) * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * intpow(cth, 2) * alpha)
	+ 0.2e1 * beta * cth * (-(a6 * a6 * h22 * g22)
	- (a6 * a6 * h21 * g21)
	- a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * intpow(sth, 2)
	+ a3 * g21 * h21 * (g22 * g22) * beta * beta * intpow(cth, 2)
	- a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * intpow(sth, 2)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a3 * (g21 * g21) * g22 * h22 * beta * beta * intpow(cth, 2)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * intpow(sth, 2)
	+ a3 * intpow(g21, 3) * h21 * beta * beta * intpow(cth, 2)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * intpow(sth, 2));

	F3 = -0.2e1 * beta * cth * (0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * h21 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g21 * h21 * beta * intpow(cth, 2) * alpha
	- 0.2e1 * a3 * g12 * h22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g11 * h21 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.2e1 * a3 * g11 * h21 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g22 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g22 * h22 * beta * intpow(cth, 2) * alpha)
	+ 0.4e1 * alpha * sth * (-0.2e1 * a6 * a6 * h22 * g22
	- 0.2e1 * a6 * a6 * h21 * g21
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	- 0.2e1 * a3 * intpow(g22, 3) * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g21 * h21 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * h21 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g21 * g21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a3 * g21 * h21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * intpow(g21, 3) * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g21 * g21 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g22 * g22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g11 * h21 * g21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g12 * h22 * g22 * g22 * alpha * alpha * sth * cth)
	+ 0.2e1 * beta * cth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * intpow(sth, 2) * beta
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * intpow(cth, 2) * alpha);

	G0 = -beta * beta * intpow(cth, 2) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G1 = -beta * beta * intpow(cth, 2) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta)
	+ 0.4e1 * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G2 = -beta * beta * intpow(cth, 2) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	+ 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta)
	- (-0.2e1 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * alpha * alpha * intpow(sth, 2)) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G3 = -beta * beta * intpow(cth, 2) * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha)
	+ 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- (-0.2e1 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * alpha * alpha * intpow(sth, 2)) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta)
	- 0.4e1 * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G4 = -0.2e1 * beta * beta * intpow(cth, 2) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth)
	+ 0.4e1 * beta * cth * alpha * sth * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha)
	- (-0.2e1 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * alpha * alpha * intpow(sth, 2)) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * intpow(sth, 2)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * intpow(cth, 2) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * intpow(cth, 2) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * intpow(sth, 2) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * intpow(cth, 2) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * intpow(sth, 2) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * intpow(sth, 2) * beta);

	H0 = -beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2));

	H1 = -beta * beta * sth * cth * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta)
	- (-0.2e1 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * alpha * intpow(cth, 2) * beta) * (a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2));

	H2 = -beta * beta * sth * cth * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * alpha * intpow(cth, 2) * beta) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2));

	H3 = -beta * beta * sth * cth * (0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta)
	- (-0.2e1 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * alpha * intpow(cth, 2) * beta) * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta)
	- (-0.2e1 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * beta * intpow(sth, 2) * alpha) * (a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2));

	H4 = -0.2e1 * beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2))
	- (-0.2e1 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * alpha * intpow(cth, 2) * beta) * (0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * beta * intpow(sth, 2) * alpha) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * intpow(cth, 2) * beta);

	J0 = -beta * cth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * intpow(sth, 2)
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * intpow(g12, 3) * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g11, 3) * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J1 = -beta * cth * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 0.2e1 * alpha * sth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * intpow(sth, 2)
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * intpow(g12, 3) * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g11, 3) * h21 * beta * beta * intpow(sth, 2)
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J2 = -beta * cth * (-(2 * a4 * a6 * g12 * h22)
	- (2 * a4 * a6 * g11 * h21)
	- 0.2e1 * a5 * g12 * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * (g12 * g12) * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * intpow(g12, 3) * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * g22 * g12 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * intpow(g11, 3) * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * (g11 * g11) * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth
	+ 0.4e1 * a5 * (g11 * g11) * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * (g11 * g11) * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * (g11 * g11) * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * (g12 * g12) * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * (g12 * g12) * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a5 * (g12 * g12) * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * alpha * intpow(sth, 2))
	+ 0.2e1 * alpha * sth * (0.4e1 * a5 * (g11 * g11) * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * (g11 * g11) * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * (g11 * g11) * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * (g11 * g11) * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * (g12 * g12) * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * (g12 * g12) * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * (g12 * g12) * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * (g12 * g12) * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * (g12 * g12) * g22 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * (g11 * g11) * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ beta * cth * (-(a4 * a6 * g12 * h22)
	- (a4 * a6 * g11 * h21)
	+ a5 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * (g12 * g12) * g11 * h21 * beta * beta * intpow(sth, 2)
	- a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * intpow(g12, 3) * h22 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g11, 3) * h21 * beta * beta * intpow(sth, 2)
	+ a5 * (g11 * g11) * g12 * h22 * beta * beta * intpow(sth, 2)
	- a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ a5 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth);

	J3 = -beta * cth * (-0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g11 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 0.2e1 * alpha * sth * (-0.2e1 * a4 * a6 * g12 * h22
	- 0.2e1 * a4 * a6 * g11 * h21
	- 0.2e1 * a5 * g12 * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g12 * g12 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h21 * beta * beta * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * intpow(g12, 3) * h22 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * g22 * g12 * h22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * intpow(g11, 3) * h21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g11 * g11 * g12 * h22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g22 * h22 * beta * beta * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * g11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g11 * g11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * intpow(sth, 2)
	+ 0.8e1 * a5 * g12 * g12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * alpha * intpow(sth, 2))
	+ beta * cth * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g12, 3) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g11, 3) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth);

	K0 = -beta * sth * (a5 * intpow(g21, 3) * h11 * beta * beta * intpow(cth, 2)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g22, 3) * h12 * beta * beta * intpow(cth, 2)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K1 = -beta * sth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * intpow(sth, 2) * alpha)
	- 0.2e1 * alpha * cth * (a5 * intpow(g21, 3) * h11 * beta * beta * intpow(cth, 2)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g22, 3) * h12 * beta * beta * intpow(cth, 2)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K2 = -beta * sth * (-0.2e1 * a5 * intpow(g21, 3) * h11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g21 * h11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * intpow(g22, 3) * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g22 * h12 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a4 * a6 * h11 * g21
	- 0.2e1 * a4 * a6 * h12 * g22
	+ 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * intpow(sth, 2))
	- 0.2e1 * alpha * cth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * intpow(sth, 2) * alpha)
	+ beta * sth * (a5 * intpow(g21, 3) * h11 * beta * beta * intpow(cth, 2)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * h11 * beta * beta * intpow(cth, 2)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * intpow(sth, 2)
	+ a5 * g22 * h12 * beta * beta * intpow(cth, 2)
	+ a5 * intpow(g22, 3) * h12 * beta * beta * intpow(cth, 2)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K3 = -beta * sth * (-0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g12 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g22 * h12 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h11 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g21 * h11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * intpow(sth, 2) * alpha)
	- 0.2e1 * alpha * cth * (-0.2e1 * a5 * intpow(g21, 3) * h11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g21 * h11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * g22 * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * intpow(g22, 3) * h12 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g22 * h12 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a4 * a6 * h11 * g21
	- 0.2e1 * a4 * a6 * h12 * g22
	+ 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * intpow(cth, 2)
	+ 0.4e1 * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * alpha * intpow(sth, 2)
	+ 0.4e1 * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * alpha * intpow(cth, 2)
	+ 0.8e1 * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * intpow(sth, 2))
	+ beta * sth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * intpow(sth, 2) * alpha);

L0 = a4 * a6 * g11 * g21
	+ a4 * a6 * g12 * g22
	+ a4 * a6
	- a5 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a5 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ a5 * g11 * g11 * g22 * g22 * beta * beta * sth * cth
	- a5 * beta * beta * intpow(sth, 2)
	- a5 * beta * beta * intpow(cth, 2)
	- a5 * intpow(g12, 3) * g22 * beta * beta * intpow(sth, 2)
	- a5 * g12 * g12 * g11 * g21 * beta * beta * intpow(sth, 2)
	- a5 * intpow(g21, 3) * g11 * beta * beta * intpow(cth, 2)
	- a5 * g22 * g12 * beta * beta * intpow(cth, 2)
	- a5 * g21 * g11 * beta * beta * intpow(sth, 2)
	- a5 * intpow(g22, 3) * g12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * sth * cth
	- a5 * g22 * g12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * beta * beta * intpow(cth, 2)
	+ a5 * g12 * g12 * beta * beta * sth * cth
	- a5 * g12 * g12 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * g11 * g11 * g21 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * sth * cth
	+ 0.2e1 * a5 * beta * beta * sth * cth
	- a5 * g22 * g12 * g21 * g21 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g12 * g12 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * g11 * beta * beta * sth * cth
	- a5 * g11 * g11 * g12 * g22 * beta * beta * intpow(sth, 2)
	- a5 * intpow(g11, 3) * g21 * beta * beta * intpow(sth, 2)
	+ a5 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g22 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * g22 * g22 * beta * beta * intpow(cth, 2)
	- a5 * g11 * g11 * beta * beta * intpow(sth, 2);

L1 = 0.4e1 * a5 * intpow(g22, 3) * g12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * intpow(g21, 3) * g11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * g22 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g21 * g11 * g12 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g21 * g11 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g21 * g11 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * g12 * g11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * intpow(g12, 3) * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * g11 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g11 * g11 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g12 * g12 * g22 * g22 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g12 * g12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g12 * alpha * intpow(cth, 2) * beta
	- 0.4e1 * a5 * g22 * g12 * beta * intpow(sth, 2) * alpha
	- 0.2e1 * a5 * g21 * g21 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * g21 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g11 * g11 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.4e1 * a5 * g22 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * g12 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * g12 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g11 * g11 * g22 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.4e1 * a5 * g21 * g11 * g22 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * intpow(g11, 3) * g21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * g21 * alpha * intpow(cth, 2) * beta
	+ 0.2e1 * a5 * g22 * g22 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g22 * g22 * beta * intpow(sth, 2) * alpha
	- 0.4e1 * a5 * g11 * g11 * g12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * g22 * beta * intpow(sth, 2) * alpha
	+ 0.2e1 * a5 * g12 * g12 * g21 * g21 * alpha * intpow(cth, 2) * beta
	- 0.2e1 * a5 * g12 * g12 * g21 * g21 * beta * intpow(sth, 2) * alpha;

L2 = (2 * a4 * a6 * g11 * g21)
	+ (2 * a4 * a6 * g12 * g22)
	+ (2 * a4 * a6)
	+ 0.2e1 * a5 * (g22 * g22) * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * (g21 * g21) * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * (g11 * g11) * (g22 * g22) * beta * beta * sth * cth
	+ 0.2e1 * a5 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * intpow(g12, 3) * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * (g12 * g12) * g11 * g21 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * intpow(g21, 3) * g11 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * intpow(g22, 3) * g12 * beta * beta * intpow(cth, 2)
	- 0.4e1 * a5 * g21 * g11 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * intpow(sth, 2)
	- 0.4e1 * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * intpow(cth, 2)
	- 0.2e1 * a5 * (g12 * g12) * beta * beta * sth * cth
	+ 0.2e1 * a5 * (g12 * g12) * beta * beta * intpow(sth, 2)
	- 0.4e1 * a5 * (g11 * g11) * (g21 * g21) * beta * beta * sth * cth
	- 0.4e1 * a5 * g22 * g12 * beta * beta * sth * cth
	- 0.4e1 * a5 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * (g21 * g21) * beta * beta * intpow(cth, 2)
	- 0.4e1 * a5 * (g12 * g12) * (g22 * g22) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g12 * g12) * (g21 * g21) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g11 * g11) * beta * beta * sth * cth
	+ 0.2e1 * a5 * (g11 * g11) * g12 * g22 * beta * beta * intpow(sth, 2)
	+ 0.2e1 * a5 * intpow(g11, 3) * g21 * beta * beta * intpow(sth, 2)
	- 0.2e1 * a5 * (g21 * g21) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g22 * g22) * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * (g22 * g22) * beta * beta * intpow(cth, 2)
	+ 0.2e1 * a5 * (g11 * g11) * beta * beta * intpow(sth, 2)
	- 0.4e1 * a5 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * intpow(g22, 3) * g12 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * (g21 * g21) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a5 * (g11 * g11) * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g22 * g12 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * (g11 * g11) * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * (g22 * g22) * alpha * alpha * intpow(sth, 2)
	- 0.8e1 * a5 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * intpow(g21, 3) * g11 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * g22 * g12 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * (g12 * g12) * g11 * g21 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * g21 * g11 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * intpow(g12, 3) * g22 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * (g11 * g11) * alpha * alpha * cth * sth
	- 0.8e1 * a5 * (g12 * g12) * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g12 * g12) * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a5 * g21 * g11 * g12 * g22 * alpha * alpha * cth * sth
	- 0.8e1 * a5 * g22 * g12 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g21 * g11 * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * (g11 * g11) * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g12 * g12) * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * intpow(g11, 3) * g21 * alpha * alpha * intpow(cth, 2)
	- 0.4e1 * a5 * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g11 * g11) * g12 * g22 * alpha * alpha * intpow(cth, 2)
	- 0.8e1 * a5 * g21 * g11 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g22 * g12 * (g21 * g21) * alpha * alpha * intpow(sth, 2)
	- 0.4e1 * a5 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g21 * g11 * (g22 * g22) * alpha * alpha * intpow(sth, 2);
}

}

#endif // !pose_from_point_tangents_2_hxx_


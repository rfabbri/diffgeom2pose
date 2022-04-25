//% This is usually called from pose_from_point_tangents_root_find_function_any.m
//%
//% Input: gama1, gama2, tgt1, tgt2, Gama1, Gama2, Tgt1, Tgt2
//% ----------------------------------------------------------------------------
//
//% dummy input ---------------------------
//%gama1 = rand(1,3);
//%gama1(3) = 1;
//%gama2 = rand(1,3);
//%gama2(3) = 1;
//%tgt1 = rand(1,3);
//%tgt1(3) = 0;
//%tgt2 = rand(1,3);
//%tgt2(3) = 0;
//
//%Gama1 = rand(1,3);
//%Gama2 = rand(1,3);
//%Tgt1 = rand(1,3);
//%Tgt1 = Tgt1/norm(Tgt1);
//%Tgt2 = rand(1,3);
//%Tgt2 = Tgt2/norm(Tgt2);
//
//%t=1;
//
//% !dummy input ---------------------------
//
//
//% Precomputed, t-independent terms:

#include "common.hxx"

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
	constexpr double PI = 3.141592653589793;

	static double g11, g12, g21, g22, h11, h12, h21, h22;
	static double V[3], buf[3];
	static double a1, a2, a3, a4, a5, a6;
	static double t1, t2, t4, t5, t6, t7, t8;
	static double t11, t14, t15, t21, t24, t25;
	static double den1, den2;
	static double sth, cth;

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

	#if 0
	A0 = a4 * a4 * g12 * g12 + a4 * a4 * g11 * g11 + a4 * a4 + 2.0 * a2 * pow(g11,
	3.0) * g21 * beta * beta * sin(theta) * cos(theta) + 2.0 * a2 * g21 * g11 *
	g12 * g12 * beta * beta * sin(theta) * cos(theta) - 2.0 * a2 * g11 * g11 * g12
	* g12 * beta * beta * pow(sin(theta), 2.0) - a2 * pow(g12, 4.0) * beta *
	beta * pow(sin(theta), 2.0) - a2 * g21 * g21 * g11 * g11 * beta * beta *
	pow(cos(theta), 2.0) + 2.0 * a2 * g12 * g12 * beta * beta * sin(theta) *
	cos(theta) + 2.0 * a2 * g11 * g11 * beta * beta * sin(theta) * cos(theta) +
	2.0 * a2 * g11 * g11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - a2
	* beta * beta * pow(cos(theta), 2.0) + 2.0 * a2 * pow(g12, 3.0) * g22 *
	beta * beta * sin(theta) * cos(theta) - a2 * pow(g11, 4.0) * beta * beta *
	pow(sin(theta), 2.0) - 2.0 * a2 * g11 * g11 * beta * beta * pow(sin(theta),
	2.0) - 2.0 * a2 * g12 * g12 * beta * beta * pow(sin(theta), 2.0) + 2.0 *
	a2 * beta * beta * sin(theta) * cos(theta) - 2.0 * a2 * g21 * g11 * g22 * g12
	* beta * beta * pow(cos(theta), 2.0) - a2 * beta * beta * pow(sin(theta),
	2.0) + 2.0 * a2 * g21 * g11 * beta * beta * sin(theta) * cos(theta) - a2 *
	g22 * g22 * g12 * g12 * beta * beta * pow(cos(theta), 2.0) - 2.0 * a2 * g22
	* g12 * beta * beta * pow(cos(theta), 2.0) - 2.0 * a2 * g21 * g11 * beta *
	beta * pow(cos(theta), 2.0) + 2.0 * a2 * g22 * g12 * beta * beta *
	sin(theta) * cos(theta);

	A1 = 0.4e1 * a2 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g21 * g11 * g11 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.8e1 * a2 * g21 * g11 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a2 * g12 * g12 * beta * pow(sin(theta), 0.2e1) *
	alpha + 0.4e1 * a2 * g22 * g12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 *
	a2 * g22 * g12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g22 * g22
	* g12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 *
	g12 * g12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g11 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a2 * g21 * g11 * g12 * g12 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.4e1 * a2 * g21 * g11 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.8e1 * a2 * g11 * g11 * g12 * g12 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a2 * pow(g11, 0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.8e1 *
	a2 * g11 * g11 * alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g21 * g11
	* g22 * g12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g12,
	0.3e1) * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * pow(g12,
	0.3e1) * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * pow(g12,
	0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.8e1 * a2 * g12 * g12 * alpha
	* sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g22 * g12 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a2 * g11 * g11 * g22 * g12 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * g22 * g12 * beta * pow(sin(theta),
	0.2e1) * alpha;

	A2 =  (2 * a4 * a4 * g12 * g12) +  (2 * a4 * a4 * g11 * g11) +
	 (2 * a4 * a4) + 0.2e1 * a2 *   pow( g12,
	4) * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 *
	pow( g11,  4) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 *
	a2 *  (g11 * g11) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 *
	 (g12 * g12) * beta * beta * pow(sin(theta), 0.2e1) - 0.4e1 * a2 * beta
	* beta * sin(theta) * cos(theta) + 0.2e1 * a2 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a2 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 *
	 g11 *  (g12 * g12) * beta * beta * sin(theta) * cos(theta) +
	0.2e1 * a2 * g21 * g21 *  (g11 * g11) * beta * beta * pow(cos(theta),
	0.2e1) - 0.4e1 * a2 *  (g12 * g12) * beta * beta * sin(theta) *
	cos(theta) + 0.4e1 * a2 * g21 *  g11 * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a2 * g22 * g22 *  (g12 * g12) * beta * beta *
	pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g22 *  g12 * beta * beta *
	sin(theta) * cos(theta) - 0.4e1 * a2 *  (g11 * g11) * beta * beta *
	sin(theta) * cos(theta) - 0.4e1 * a2 * g21 *  g11 * beta * beta *
	sin(theta) * cos(theta) + 0.4e1 * a2 *  (g11 * g11) *  (g12 *
	g12) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 *  g12 *
	beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 *   pow(
	g11,  3) * g21 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a2 *
	 (g11 * g11) * g22 *  g12 * beta * beta * sin(theta) *
	cos(theta) + 0.4e1 * a2 * g21 *  g11 * g22 *  g12 * beta * beta
	* pow(cos(theta), 0.2e1) - 0.4e1 * a2 *   pow( g12,
	 3) * g22 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a2 *
	  pow( g11,  4) * alpha * alpha * pow(cos(theta),
	0.2e1) - 0.8e1 * a2 *  (g11 * g11) * alpha * alpha * pow(cos(theta),
	0.2e1) - 0.4e1 * a2 *   pow( g12,  4) * alpha *
	alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a2 *  (g12 * g12) * alpha *
	alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a2 * alpha * alpha * cos(theta) *
	sin(theta) - 0.4e1 * a2 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a2 *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 * g22 *  g12 * alpha
	* alpha * cos(theta) * sin(theta) - 0.4e1 * a2 * g21 * g21 *  (g11 *
	g11) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 *  (g12 *
	g12) * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a2 * g21 *  g11
	* alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a2 * g22 * g22 *
	(g12 * g12) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 * g21 *
	 g11 * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a2 *
	(g11 * g11) *  (g12 * g12) * alpha * alpha * pow(cos(theta), 0.2e1) -
	0.8e1 * a2 *  (g11 * g11) * alpha * alpha * cos(theta) * sin(theta) -
	0.8e1 * a2 * g21 *  g11 *  (g12 * g12) * alpha * alpha *
	cos(theta) * sin(theta) - 0.8e1 * a2 * g21 *  g11 * g22 *  g12 *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 *
	pow( g12,  3) * g22 * alpha * alpha * cos(theta) * sin(theta) -
	0.8e1 * a2 * g22 *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1
	* a2 *   pow( g11,  3) * g21 * alpha * alpha *
	cos(theta) * sin(theta) - 0.8e1 * a2 *  (g11 * g11) * g22 *  g12
	* alpha * alpha * cos(theta) * sin(theta);

	B0 = -0.2e1 * beta * sin(theta) * (a2 * g21 * g11 * g22 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) + a2 * pow(g12, 0.3e1) * h12 * beta * beta *
	pow(sin(theta), 0.2e1) + a2 * g21 * h11 * g22 * g12 * beta * beta *
	pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 * beta * beta *
	pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta) * cos(theta)
	- a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) +
	a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g22 *
	h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 * g12 * beta * beta * sin(theta) *
	cos(theta) - a2 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 *
	h11 * g12 * g12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * h11 * beta
	* beta * sin(theta) * cos(theta) - a2 * g11 * g11 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta *
	sin(theta) * cos(theta) - a2 * g22 * h12 * beta * beta * sin(theta) * cos(theta)
	+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 *
	g21 * h11 * g11 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta *
	beta * pow(cos(theta), 0.2e1) - a2 * g21 * g11 * g12 * h12 * beta * beta *
	sin(theta) * cos(theta) + a2 * g11 * g11 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

	B1 = -0.2e1 * beta * sin(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sin(theta) * beta
	* cos(theta) - 0.2e1 * a2 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a2 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 *
	g22 * h12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta)
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a2 * g11 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 *
	g21 * g11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 *
	g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11
	* g11 * alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * alpha * cos(theta) *
	(a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 *
	pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * h11 *
	g22 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 *
	beta * beta * pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta)
	* cos(theta) - a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta *
	pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * g12 * beta * beta *
	pow(cos(theta), 0.2e1) + a2 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) +
	a2 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 *
	g12 * beta * beta * sin(theta) * cos(theta) - a2 * g12 * h12 * beta * beta *
	sin(theta) * cos(theta) - a2 * g21 * h11 * g12 * g12 * beta * beta * sin(theta)
	* cos(theta) - a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - a2 * g11
	* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 *
	g12 * h12 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11
	* h11 * g21 * beta * beta * sin(theta) * cos(theta) - a2 * g22 * h12 * beta *
	beta * sin(theta) * cos(theta) + a2 * pow(g11, 0.3e1) * h11 * beta * beta *
	pow(sin(theta), 0.2e1) + a2 * g21 * g21 * h11 * g11 * beta * beta *
	pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) -
	a2 * g21 * g11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g11 *
	g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

	B2 = -0.2e1 * beta * sin(theta) * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha *
	alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha *
	alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g11 * h11 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * alpha * alpha * cos(theta) *
	sin(theta) + 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 * h12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) +
	0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) -
	0.2e1 * a4 * a4 * h11 * g11 - 0.2e1 * a4 * a4 * h12 * g12 - 0.2e1 * a2 * g22 *
	g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g12 * h12
	* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 *
	a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 *
	g11 * g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11
	* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g11 *
	h11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * h12
	* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sin(theta),
	0.2e1) - 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta),
	0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1
	* a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2
	* g21 * h11 * g22 * g12 * beta * beta * pow(cos(theta), 0.2e1)) - 0.4e1 * alpha
	* cos(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sin(theta) * beta
	* cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 *
	g22 * g22 * h12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 *
	g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g22 * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sin(theta) *
	beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 *
	a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * g11
	* g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 *
	g12 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 * g22 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha *
	sin(theta) * beta * cos(theta)) + 0.2e1 * beta * sin(theta) * (a2 * g21 * g11 *
	g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * pow(g12, 0.3e1) * h12 *
	beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * h11 * g22 * g12 * beta * beta
	* pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 * beta * beta *
	pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta) * cos(theta)
	- a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) +
	a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g22 *
	h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 * g12 * beta * beta * sin(theta) *
	cos(theta) - a2 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 *
	h11 * g12 * g12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * h11 * beta
	* beta * sin(theta) * cos(theta) - a2 * g11 * g11 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta *
	sin(theta) * cos(theta) - a2 * g22 * h12 * beta * beta * sin(theta) * cos(theta)
	+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 *
	g21 * h11 * g11 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta *
	beta * pow(cos(theta), 0.2e1) - a2 * g21 * g11 * g12 * h12 * beta * beta *
	sin(theta) * cos(theta) + a2 * g11 * g11 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

	B3 = -0.2e1 * beta * sin(theta) * (-0.2e1 * a2 * g11 * g11 * g22 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha *
	sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sin(theta) * beta
	* cos(theta) + 0.2e1 * a2 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta -
	0.2e1 * a2 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 *
	g22 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * h11 * g22 * g12 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g12 * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha *
	sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g21 * h11 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta)
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta -
	0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a2 * g11 * h11 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 *
	g21 * g11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 *
	g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g21 * h11
	* g11 * alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * alpha * cos(theta) *
	(0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a2 * g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 *
	g11 * g11 * g22 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 *
	g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 * g11 * g11 *
	h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * pow(g12,
	0.3e1) * h12 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g22 * g12 *
	alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * alpha * alpha
	* cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 * h12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) +
	0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) -
	0.2e1 * a4 * a4 * h11 * g11 - 0.2e1 * a4 * a4 * h12 * g12 - 0.2e1 * a2 * g22 *
	g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g12 * h12
	* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 *
	a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 *
	g11 * g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11
	* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g11 *
	h11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * h12
	* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sin(theta),
	0.2e1) - 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta),
	0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1
	* a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2
	* g21 * h11 * g22 * g12 * beta * beta * pow(cos(theta), 0.2e1)) + 0.2e1 * beta *
	sin(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sin(theta), 0.2e1)
	* alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) *
	beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sin(theta), 0.2e1) *
	alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sin(theta), 0.2e1) *
	alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sin(theta) * beta *
	cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 *
	g22 * g22 * h12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 *
	g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g22 * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sin(theta) *
	beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 *
	a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * g11
	* g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 *
	g12 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 * g22 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha *
	sin(theta) * beta * cos(theta));

	C0 = -beta * beta * pow(sin(theta), 0.2e1) * (-a4 * a4 * h12 * h12 + 0.2e1 * a2
	* g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 *
	g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4 * a4 * h11 *
	h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta) * cos(theta)
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2
	* g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g12 * g12
	* h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12
	* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 *
	beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta
	* beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta *
	pow(sin(theta), 0.2e1));

	C1 = -beta * beta * pow(sin(theta), 0.2e1) * (0.8e1 * a2 * g11 * h11 * g12 *
	h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 *
	h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 *
	h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 *
	h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22
	* alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 *
	alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * beta * sin(theta) * alpha *
	cos(theta) * (-a4 * a4 * h12 * h12 + 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta *
	beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta
	* sin(theta) * cos(theta) - a4 * a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 *
	g22 * beta * beta * sin(theta) * cos(theta) + a2 * g21 * g21 * h11 * h11 * beta
	* beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta *
	beta * sin(theta) * cos(theta) + a2 * g12 * g12 * h12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta *
	pow(sin(theta), 0.2e1));

	C2 = -beta * beta * pow(sin(theta), 0.2e1) * (-0.4e1 * a2 * g11 * h11 * g12 *
	h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22
	* beta * beta * sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	beta * beta * sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2
	* g21 * g21 * h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 *
	g21 * h11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 *
	g12 * g12 * h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 *
	g22 * g22 * h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21
	* h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 *
	h11 * g22 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 *
	h11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11
	* g12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 *
	g22 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 *
	h11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11
	* h11 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 +
	0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) +
	0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) +
	0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) -
	0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) -
	0.4e1 * beta * sin(theta) * alpha * cos(theta) * (0.8e1 * a2 * g11 * h11 * g12 *
	h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 *
	h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 *
	h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 *
	h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22
	* alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 *
	alpha * sin(theta) * beta * cos(theta)) - (-0.2e1 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * alpha * alpha * pow(cos(theta), 0.2e1)) * (-a4
	* a4 * h12 * h12 + 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta *
	sin(theta) * cos(theta) - a4 * a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 *
	g22 * beta * beta * sin(theta) * cos(theta) + a2 * g21 * g21 * h11 * h11 * beta
	* beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta *
	beta * sin(theta) * cos(theta) + a2 * g12 * g12 * h12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta *
	sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta *
	pow(sin(theta), 0.2e1));

	C3 = -beta * beta * pow(sin(theta), 0.2e1) * (0.4e1 * a2 * g21 * h11 * g12 * h12
	* alpha * pow(cos(theta), 0.2e1) * beta - 0.8e1 * a2 * g11 * h11 * g12 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g21 * h11 * g22 * h12 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g12 * h12 * h12 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * h12 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * h12 * g22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * g21 * h11 * h11 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g21 * h11 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * g11 * h11 * h11 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g22 * g22 * h12 * h12 *
	alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * beta * sin(theta) * alpha *
	cos(theta) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta *
	sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2 * g21 * g21 *
	h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 *
	g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g12 * g12 *
	h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g22 * g22 *
	h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 * h11 * g22
	* h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 * h11 * g22 *
	h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 *
	h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11 * g12 * h12
	* alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * g22 * h12 *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 * h11 * h11 *
	beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * h11 *
	alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 + 0.4e1 *
	a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2
	* g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 *
	g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 *
	g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12
	* g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g21 *
	g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) - (-0.2e1 * beta * beta
	* pow(sin(theta), 0.2e1) + 0.4e1 * alpha * alpha * pow(cos(theta), 0.2e1)) *
	(0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta -
	0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cos(theta), 0.2e1) * beta -
	0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sin(theta) * beta * cos(theta)) +
	0.4e1 * beta * sin(theta) * alpha * cos(theta) * (-a4 * a4 * h12 * h12 + 0.2e1 *
	a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 *
	g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4 * a4 * h11 *
	h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta) * cos(theta)
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2
	* g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g12 * g12
	* h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12
	* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 *
	beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta
	* beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta *
	pow(sin(theta), 0.2e1));

	C4 = -0.2e1 * beta * beta * pow(sin(theta), 0.2e1) * (-a4 * a4 * h12 * h12 +
	0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4
	* a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta)
	* cos(theta) + a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) +
	a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 *
	g22 * h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11
	* g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 *
	g22 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 *
	beta * beta * pow(sin(theta), 0.2e1)) - 0.4e1 * beta * sin(theta) * alpha *
	cos(theta) * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sin(theta) * beta *
	cos(theta) + 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sin(theta), 0.2e1)
	* alpha - 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1)
	* beta + 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sin(theta), 0.2e1)
	* alpha + 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) *
	beta - 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cos(theta), 0.2e1)
	* beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sin(theta), 0.2e1) *
	alpha + 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sin(theta) * beta *
	cos(theta)) - (-0.2e1 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * alpha *
	alpha * pow(cos(theta), 0.2e1)) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta *
	beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta
	* sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2 * g21 * g21 *
	h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 *
	g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g12 * g12 *
	h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g22 * g22 *
	h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 * h11 * g22
	* h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 * h11 * g22 *
	h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 *
	h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11 * g12 * h12
	* alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * g22 * h12 *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 * h11 * h11 *
	beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * h11 *
	alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 + 0.4e1 *
	a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2
	* g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 *
	g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 *
	g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12
	* g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g21 *
	g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) + 0.4e1 * beta *
	sin(theta) * alpha * cos(theta) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha *
	sin(theta) * beta * cos(theta));

	E0 = 0.2e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cos(theta) *
	sin(theta) + 0.2e1 * a3 * g12 * g22 * beta * beta * cos(theta) * sin(theta) +
	0.2e1 * a3 * g12 * pow(g22, 0.3e1) * beta * beta * cos(theta) * sin(theta) - a3
	* beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * g21 * beta * beta * pow(sin(theta),
	0.2e1) - a3 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 * g11 * g21 *
	beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g22 * g22 * beta * beta *
	cos(theta) * sin(theta) - a3 * pow(g21, 0.4e1) * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a3 * g21 * g21 * beta * beta * cos(theta) * sin(theta) - a3 *
	g12 * g12 * g22 * g22 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 * g11 *
	g21 * g21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g21 * g21 * g22
	* g22 * beta * beta * pow(cos(theta), 0.2e1) + a6 * a6 * g21 * g21 + a6 * a6 *
	g22 * g22 + 0.2e1 * a3 * g11 * pow(g21, 0.3e1) * beta * beta * cos(theta) *
	sin(theta) + 0.2e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cos(theta) *
	sin(theta) - 0.2e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + a6 * a6 - 0.2e1 * a3 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1)
	- 0.2e1 * a3 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * pow(g22,
	0.4e1) * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 * beta * beta *
	cos(theta) * sin(theta);

	E1 = -0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * sin(theta) * beta *
	cos(theta) + 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a3 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha -
	0.4e1 * a3 * g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.8e1 * a3 *
	g11 * g21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * g22 *
	alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g22 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * g21 * g12 * g22 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.8e1 * a3 * g21 * g21 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a3 * g11 * g21 * alpha * pow(sin(theta), 0.2e1) * beta +
	0.4e1 * a3 * g11 * g21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 *
	g12 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 *
	g12 * g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g21 * g22
	* g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g11 * g21 * g12 *
	g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * pow(g21,
	0.3e1) * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * pow(g21,
	0.3e1) * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * g21 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.8e1 * a3 * g22 * g22 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a3 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 *
	a3 * g21 * g21 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * pow(g22,
	0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * g21 * g22 *
	g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * beta *
	pow(cos(theta), 0.2e1) * alpha;

	E2 = -0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cos(theta) * sin(theta)
	- 0.4e1 * a3 * g12 * g22 * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 *
	g12 * pow(g22, 0.3e1) * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 *
	beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g21 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a3 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a3 * g11 *
	g21 * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 * g22 * g22 * beta *
	beta * cos(theta) * sin(theta) + 0.2e1 * a3 * pow(g21, 0.4e1) * beta * beta *
	pow(cos(theta), 0.2e1) - 0.4e1 * a3 * g21 * g21 * beta * beta * cos(theta) *
	sin(theta) + 0.2e1 * a3 * g12 * g12 * g22 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a3 * g11 * g11 * g21 * g21 * beta * beta * pow(sin(theta),
	0.2e1) + 0.4e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a6 * a6 * g21 * g21 + 0.2e1 * a6 * a6 * g22 * g22 - 0.4e1 * a3
	* g11 * pow(g21, 0.3e1) * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 *
	g11 * g21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) + 0.4e1 * a3 * g11
	* g21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a6 * a6 -
	0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) -
	0.8e1 * a3 * g12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a3 *
	g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a3 * g12 * g22 *
	alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 * g12 * pow(g22, 0.3e1) *
	alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 * g21 * g21 * g22 * g22 *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a3 * g22 * g22 * alpha * alpha
	* sin(theta) * cos(theta) - 0.8e1 * a3 * g21 * g21 * g12 * g22 * alpha * alpha *
	sin(theta) * cos(theta) - 0.8e1 * a3 * g11 * pow(g21, 0.3e1) * alpha * alpha *
	sin(theta) * cos(theta) - 0.4e1 * a3 * alpha * alpha * pow(sin(theta), 0.2e1) -
	0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) -
	0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) -
	0.8e1 * a3 * g21 * g21 * alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 *
	g11 * g21 * alpha * alpha * sin(theta) * cos(theta) - 0.4e1 * a3 * alpha * alpha
	* pow(cos(theta), 0.2e1) - 0.8e1 * a3 * g11 * g21 * g22 * g22 * alpha * alpha *
	sin(theta) * cos(theta) + 0.4e1 * a3 * g21 * g21 * beta * beta * pow(cos(theta),
	0.2e1) + 0.4e1 * a3 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 *
	a3 * pow(g22, 0.4e1) * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a3 * beta
	* beta * cos(theta) * sin(theta) - 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * alpha
	* pow(sin(theta), 0.2e1) - 0.8e1 * a3 * g21 * g21 * alpha * alpha *
	pow(sin(theta), 0.2e1) - 0.8e1 * a3 * alpha * alpha * sin(theta) * cos(theta) -
	0.4e1 * a3 * pow(g22, 0.4e1) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 *
	a3 * g22 * g22 * alpha * alpha * pow(sin(theta), 0.2e1);

	F0 = -0.2e1 * beta * cos(theta) * (-a6 * a6 * h22 * g22 - a6 * a6 * h21 * g21 -
	a3 * g11 * h21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) + a3 * g11 *
	h21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21
	* g21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * pow(g22, 0.3e1) * h22
	* beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - a3 * g21 * h21 * g12 * g22 * beta * beta * cos(theta) *
	sin(theta) + a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g12 *
	g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g21 * h21 * g22 *
	g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * g21 * g21 * g12 * h22 * beta *
	beta * cos(theta) * sin(theta) - a3 * g11 * h21 * beta * beta * cos(theta) *
	sin(theta) + a3 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 *
	g21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) - a3 * g22 * h22 * beta
	* beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta *
	beta * cos(theta) * sin(theta) + a3 * g21 * h21 * beta * beta * pow(cos(theta),
	0.2e1) + a3 * g21 * g21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - a3
	* g21 * h21 * beta * beta * cos(theta) * sin(theta) - a3 * g12 * h22 * beta *
	beta * cos(theta) * sin(theta) + a3 * g11 * g11 * h21 * g21 * beta * beta *
	pow(sin(theta), 0.2e1) + a3 * pow(g21, 0.3e1) * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 * g11 * g21 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1));

	F1 = -0.2e1 * beta * cos(theta) * (-0.4e1 * a3 * g21 * h21 * alpha * sin(theta)
	* beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sin(theta) *
	beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * alpha * pow(sin(theta), 0.2e1) *
	beta - 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 *
	a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g12
	* h22 * g22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * h22 *
	g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * g22
	* g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g21 * h21 * g22 *
	g22 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * g21 * g12 *
	h22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * g21 * g12 * h22
	* beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11 * h21 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a3 * g22 * h22 * alpha * sin(theta) * beta *
	cos(theta) + 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sin(theta), 0.2e1)
	* beta - 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cos(theta), 0.2e1) *
	alpha + 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sin(theta) * beta *
	cos(theta) + 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1)
	* beta - 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) *
	alpha + 0.4e1 * a3 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1
	* a3 * g11 * h21 * g12 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 *
	a3 * g11 * h21 * g21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3
	* pow(g22, 0.3e1) * h22 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a3 *
	g11 * h21 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11
	* h21 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 * h22 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha *
	sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g22 * h22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta),
	0.2e1) * alpha) + 0.4e1 * alpha * sin(theta) * (-a6 * a6 * h22 * g22 - a6 * a6 *
	h21 * g21 - a3 * g11 * h21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) +
	a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 *
	g11 * h21 * g21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * pow(g22,
	0.3e1) * h22 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * h22 * beta *
	beta * pow(cos(theta), 0.2e1) - a3 * g21 * h21 * g12 * g22 * beta * beta *
	cos(theta) * sin(theta) + a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g21 *
	h21 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * g21 * g21 * g12 *
	h22 * beta * beta * cos(theta) * sin(theta) - a3 * g11 * h21 * beta * beta *
	cos(theta) * sin(theta) + a3 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) - a3 * g22
	* h22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 * h22 * g22 *
	g22 * beta * beta * cos(theta) * sin(theta) + a3 * g21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 * g21 * g21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - a3 * g21 * h21 * beta * beta * cos(theta) * sin(theta)
	- a3 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + a3 * g11 * g11 * h21
	* g21 * beta * beta * pow(sin(theta), 0.2e1) + a3 * pow(g21, 0.3e1) * h21 * beta
	* beta * pow(cos(theta), 0.2e1) + a3 * g11 * g21 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1));

	F2 = -0.2e1 * beta * cos(theta) * (- (2 * a6 * a6 * h22 * g22) -  (2 * a6 * a6 *
	h21 * g21) + 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * beta * beta * cos(theta) *
	sin(theta) - 0.2e1 * a3 * g11 *  h21 * g12 *  g22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 *   pow( g22,  3) *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 *  g22 *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * g12 *  h22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g12 * g12 *  h22 *  g22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 *  g21 *  h21 *  (g22 * g22) * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  (g21 * g21) * g12 *  h22 * beta * beta *
	cos(theta) * sin(theta) + 0.2e1 * a3 * g11 *  h21 * beta * beta * cos(theta) *
	sin(theta) - 0.2e1 * a3 * g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) +
	0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) +
	0.2e1 * a3 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) + 0.4e1 * a3 *
	g12 *  h22 *  (g22 * g22) * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 *
	g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 *  (g21 * g21) *
	g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  g21 *  h21 *
	beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g12 *  h22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * g11 * g11 *  h21 *  g21 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 *   pow( g21,  3) *  h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 *  g21 * g12 *  h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  (g21 * g21) *  g22 *  h22 * alpha * alpha
	* pow(sin(theta), 0.2e1) + 0.4e1 * a3 *   pow( g21,  3) *  h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 * alpha * alpha * sin(theta) *
	cos(theta) + 0.4e1 * a3 * g12 *  h22 * alpha * alpha * sin(theta) * cos(theta) +
	0.4e1 * a3 * g12 * g12 *  h22 *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a3 *  (g21 * g21) * g12 *  h22 * alpha * alpha * sin(theta) * cos(theta)
	+ 0.4e1 * a3 * g11 *  h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3
	* g11 *  h21 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 *  g21 *
	h21 * g12 *  g22 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g11 *
	g21 * g12 *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 *
	h21 * g12 *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 *
	h21 *  (g22 * g22) * alpha * alpha * sin(theta) * cos(theta) + 0.8e1 * a3 * g11
	*  h21 *  (g21 * g21) * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 *
	g22 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 *
	(g22 * g22) * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *   pow( g22,
	3) *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g11 *
	h21 *  g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 *  h22 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 *  g21 *  g22 *  h22 *
	alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 *  g22 *  h22 * alpha *
	alpha * sin(theta) * cos(theta) + 0.8e1 * a3 * g12 *  h22 *  (g22 * g22) * alpha
	* alpha * sin(theta) * cos(theta)) + 0.4e1 * alpha * sin(theta) * (-0.4e1 * a3 *
	g21 *  h21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 *  (g21 * g21)
	*  g22 *  h22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 *   pow(
	g21,  3) *  h21 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  g21 *
	h21 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  g21 *  h21 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g12 *  h22 * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g12 *  h22 *  g22 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 *  h22 *  (g22 * g22) * alpha
	* pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 *  h22 *  (g22 * g22) * beta
	* pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 *  g21 *  h21 *  (g22 * g22) *
	alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  (g21 * g21) * g12 *  h22
	* alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  (g21 * g21) * g12 *  h22
	* beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11 *  h21 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 *  h21 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a3 *  g22 *  h22 * alpha * sin(theta) * beta *
	cos(theta) + 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * alpha * pow(sin(theta),
	0.2e1) * beta - 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * beta * pow(cos(theta),
	0.2e1) * alpha + 0.4e1 * a3 * g11 *  g21 * g12 *  h22 * alpha * sin(theta) *
	beta * cos(theta) + 0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 *  h22 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a3 * g11 *  h21 * g12 *  g22 * alpha * sin(theta) *
	beta * cos(theta) - 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 *   pow( g22,  3) *  h22 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 *  h21 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 *  h22 * beta * pow(cos(theta),
	0.2e1) * alpha + 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g11 *  h21 *  g21 * alpha *
	sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  g22 *  h22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  g22 *  h22 * beta *
	pow(cos(theta), 0.2e1) * alpha) + 0.2e1 * beta * cos(theta) * (- (a6 * a6 * h22
	* g22) -  (a6 * a6 * h21 * g21) - a3 * g11 *  h21 *  (g22 * g22) * beta * beta *
	cos(theta) * sin(theta) + a3 * g11 *  h21 * g12 *  g22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 *  h21 *  (g21 * g21) * beta * beta *
	cos(theta) * sin(theta) + a3 *   pow( g22,  3) *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 *  g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1)
	- a3 *  g21 *  h21 * g12 *  g22 * beta * beta * cos(theta) * sin(theta) + a3 *
	g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g12 * g12 *  h22 *  g22
	* beta * beta * pow(sin(theta), 0.2e1) + a3 *  g21 *  h21 *  (g22 * g22) * beta
	* beta * pow(cos(theta), 0.2e1) - a3 *  (g21 * g21) * g12 *  h22 * beta * beta *
	cos(theta) * sin(theta) - a3 * g11 *  h21 * beta * beta * cos(theta) *
	sin(theta) + a3 * g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 *
	g21 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) - a3 *  g22 *  h22 *
	beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 *  h22 *  (g22 * g22) *
	beta * beta * cos(theta) * sin(theta) + a3 *  g21 *  h21 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 *  (g21 * g21) *  g22 *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) - a3 *  g21 *  h21 * beta * beta * cos(theta) *
	sin(theta) - a3 * g12 *  h22 * beta * beta * cos(theta) * sin(theta) + a3 * g11
	* g11 *  h21 *  g21 * beta * beta * pow(sin(theta), 0.2e1) + a3 *   pow( g21,
	3) *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g11 *  g21 * g12 *  h22
	* beta * beta * pow(sin(theta), 0.2e1));

	F3 = -0.2e1 * beta * cos(theta) * (0.4e1 * a3 * g21 * h21 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sin(theta) *
	beta * cos(theta) - 0.2e1 * a3 * g21 * h21 * alpha * pow(sin(theta), 0.2e1) *
	beta + 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 *
	a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * g12
	* h22 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * h22 *
	g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * h22 * g22
	* g22 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * h21 * g22 *
	g22 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g21 * g21 * g12 *
	h22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g21 * g21 * g12 * h22
	* beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 * h21 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * alpha * sin(theta) *
	beta * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * sin(theta) * beta *
	cos(theta) - 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sin(theta), 0.2e1)
	* beta + 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cos(theta), 0.2e1) *
	alpha - 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sin(theta) * beta *
	cos(theta) - 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1)
	* beta + 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) *
	alpha - 0.4e1 * a3 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1
	* a3 * g11 * h21 * g12 * g22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 *
	a3 * g11 * h21 * g21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3
	* pow(g22, 0.3e1) * h22 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 *
	g11 * h21 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11
	* h21 * alpha * pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g12 * h22 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g22 * h22 * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta),
	0.2e1) * alpha) + 0.4e1 * alpha * sin(theta) * (-0.2e1 * a6 * a6 * h22 * g22 -
	0.2e1 * a6 * a6 * h21 * g21 + 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * pow(g22, 0.3e1) * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * h22 * beta * beta * pow(cos(theta),
	0.2e1) + 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cos(theta) *
	sin(theta) - 0.2e1 * a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a3 * g21 * h21 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) +
	0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) +
	0.2e1 * a3 * g11 * h21 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 *
	g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * g21 * g22
	* h22 * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g22 * h22 * beta *
	beta * cos(theta) * sin(theta) + 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta *
	beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g21 * g21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a3 * g21 * h21 * beta * beta * cos(theta) *
	sin(theta) + 0.2e1 * a3 * g12 * h22 * beta * beta * cos(theta) * sin(theta) -
	0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a3 * pow(g21, 0.3e1) * h21 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.4e1 * a3 * g21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *
	g21 * g21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *
	pow(g21, 0.3e1) * h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *
	g21 * h21 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g12 * h22 *
	alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * g22 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g21 * g21 * g12 * h22 *
	alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g11 * h21 * alpha * alpha
	* pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * alpha * alpha * sin(theta) *
	cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sin(theta) *
	cos(theta) + 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g22 * g22 * alpha * alpha *
	sin(theta) * cos(theta) + 0.8e1 * a3 * g11 * h21 * g21 * g21 * alpha * alpha *
	sin(theta) * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha *
	sin(theta) * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * alpha * sin(theta) *
	cos(theta) + 0.8e1 * a3 * g12 * h22 * g22 * g22 * alpha * alpha * sin(theta) *
	cos(theta)) + 0.2e1 * beta * cos(theta) * (-0.4e1 * a3 * g21 * h21 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha *
	sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta),
	0.2e1) * alpha + 0.2e1 * a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta -
	0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha -
	0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sin(theta) * beta * cos(theta) +
	0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta -
	0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) * alpha -
	0.2e1 * a3 * g11 * h21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 *
	g11 * h21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * g12 * g22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * h21 * g12 * g22 *
	beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g21 * g12 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g11 * g21 * g22 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 * g21 * g22 * h22 *
	beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha *
	sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 * h21 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 * h22 * beta * pow(cos(theta),
	0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * pow(sin(theta),
	0.2e1) * beta + 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * pow(sin(theta),
	0.2e1) * beta + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sin(theta) * beta *
	cos(theta) + 0.2e1 * a3 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) * beta -
	0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * alpha);



	G0 = -beta * beta * pow(cos(theta), 0.2e1) * (-a6 * a6 * h21 * h21 - a6 * a6 *
	h22 * h22 + a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 *
	g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11
	* h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * g21 * g21 * h21
	* h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 *
	g22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22
	* beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 *
	beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta
	* pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta *
	cos(theta) * sin(theta));

	G1= -beta * beta * pow(cos(theta), 0.2e1) * (-0.4e1 * a3 * g21 * h21 * g12 * h22
	* beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 * g21 *
	beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 *
	alpha * pow(sin(theta), 0.2e1) * beta) + 0.4e1 * beta * cos(theta) * alpha *
	sin(theta) * (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12 * g12 * h22
	* h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 *
	h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 * h21 * beta *
	beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta
	* cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta *
	cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta *
	cos(theta) * sin(theta));

	G2 = -beta * beta * pow(cos(theta), 0.2e1) * (-0.4e1 * a3 * g11 * h21 * g12 *
	h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22
	* alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 *
	beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta
	* beta * cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 *
	g21 * h21 * g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 *
	g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 *
	h22 * h22 + 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) *
	cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) *
	sin(theta) + 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) - 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha *
	sin(theta) * cos(theta) - 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta *
	cos(theta) * sin(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta *
	cos(theta) * sin(theta) + 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha *
	sin(theta) * cos(theta)) + 0.4e1 * beta * cos(theta) * alpha * sin(theta) *
	(-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) * alpha -
	0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha +
	0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) * beta -
	0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha +
	0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta -
	0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * alpha -
	0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sin(theta) * beta * cos(theta) -
	0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sin(theta) * beta * cos(theta) +
	0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) +
	0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta +
	0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sin(theta), 0.2e1) * beta) -
	(-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha *
	pow(sin(theta), 0.2e1)) * (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12
	* g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 *
	h21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 *
	h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21
	* beta * beta * cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta *
	beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta
	* cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta *
	cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta *
	cos(theta) * sin(theta));

	G3 = -beta * beta * pow(cos(theta), 0.2e1) * (0.4e1 * a3 * g22 * g22 * h22 * h22
	* alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * g21 * h21 * h21 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 *
	beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * g11 * h21 * h21 *
	alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * g12 * h22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.8e1 * a3 * g11 * h21 * g12 * h22 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * h21 * h21 * g21 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g21 * h21 * g12 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g21 * h21 * g12 * h22 *
	beta * pow(cos(theta), 0.2e1) * alpha) + 0.4e1 * beta * cos(theta) * alpha *
	sin(theta) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 * g21 * h21 *
	g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 * g11 * g11 *
	h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 * h22 * h22 +
	0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) * cos(theta) +
	0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) +
	0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) -
	0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sin(theta) * cos(theta) -
	0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) +
	0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) +
	0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) +
	0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sin(theta) * cos(theta)) -
	(-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha *
	pow(sin(theta), 0.2e1)) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta *
	pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha *
	pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta *
	pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha *
	sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha *
	sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha *
	sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha *
	pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha *
	pow(sin(theta), 0.2e1) * beta) - 0.4e1 * beta * cos(theta) * alpha * sin(theta)
	* (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12 * g12 * h22 * h22 *
	beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta
	* beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta *
	cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta *
	cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta *
	cos(theta) * sin(theta));

	G4 = -0.2e1 * beta * beta * pow(cos(theta), 0.2e1) * (-a6 * a6 * h21 * h21 - a6
	* a6 * h22 * h22 + a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sin(theta),
	0.2e1) + a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + a3
	* g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 *
	g12 * h22 * h22 * g22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21
	* h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g21 *
	h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 *
	h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22
	* beta * beta * cos(theta) * sin(theta)) + 0.4e1 * beta * cos(theta) * alpha *
	sin(theta) * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cos(theta), 0.2e1)
	* alpha - 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sin(theta), 0.2e1) *
	beta + 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) *
	alpha - 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sin(theta) * beta *
	cos(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sin(theta) * beta *
	cos(theta) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cos(theta), 0.2e1)
	* alpha - 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) *
	beta - 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sin(theta) * beta *
	cos(theta) - 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sin(theta), 0.2e1)
	* beta - 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) *
	beta + 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) *
	alpha) - (-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha
	* pow(sin(theta), 0.2e1)) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta *
	cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 * g21 * h21 *
	g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 * g11 * g11 *
	h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 * h22 * h22 +
	0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) * cos(theta) +
	0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) +
	0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) -
	0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sin(theta) * cos(theta) -
	0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) +
	0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) +
	0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) +
	0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sin(theta) * cos(theta)) -
	0.4e1 * beta * cos(theta) * alpha * sin(theta) * (-0.4e1 * a3 * g21 * h21 * g12
	* h22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 *
	g21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21
	* alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 *
	alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 *
	beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 *
	alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 *
	alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 *
	alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 *
	alpha * pow(sin(theta), 0.2e1) * beta);

	H0 = -beta * beta * sin(theta) * cos(theta) * (a5 * g11 * g11 * h11 * h21 * beta
	* beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 * h12 * h22 - a5
	* g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g22 * h12
	* g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 * g22 * h12 * h22
	* beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11 * g21 * h21 *
	beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22 * beta * beta
	* pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta * beta * sin(theta)
	* cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cos(theta), 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g12
	* h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * g22
	* h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g12 * h12 * g22 *
	h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 * g21 * h21 * beta
	* beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1));

	H1 = -beta * beta * sin(theta) * cos(theta) * (0.4e1 * a5 * g11 * g11 * h11 *
	h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 *
	h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 *
	h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21
	* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * pow(sin(theta), 0.2e1)
	* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (a5 * g11 * g11 * h11
	* h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 *
	h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) +
	a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 *
	g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11
	* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22
	* beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta *
	beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) +
	a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *
	g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 *
	g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 *
	beta * beta * pow(sin(theta), 0.2e1));

	H2 = -beta * beta * sin(theta) * cos(theta) * (-0.2e1 * a5 * g11 * g11 * h11 *
	h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 * h11 * h21 - 0.2e1
	* a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha *
	cos(theta) * sin(theta)) - (-0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (0.4e1 * a5 * g11 * g11 * h11 *
	h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 *
	h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 *
	h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21
	* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * beta * sin(theta) *
	cos(theta) - 0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (a5 * g11 * g11
	* h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 *
	a6 * h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) *
	cos(theta) + a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) +
	a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *
	g11 * h11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 *
	g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 *
	beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta
	* pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta)
	* cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5
	* g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12
	* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22
	* beta * beta * pow(sin(theta), 0.2e1));

	H3 = -beta * beta * sin(theta) * cos(theta) * (0.4e1 * a5 * g21 * g21 * h11 *
	h21 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * h11 * g12 *
	h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 *
	h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * h12 * g11 * h21
	* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * pow(sin(theta), 0.2e1)
	* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (-0.2e1 * a5 * g11 *
	g11 * h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 * h11 *
	h21 - 0.2e1 * a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta *
	beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta *
	beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta
	* pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha *
	cos(theta) * sin(theta)) - (-0.2e1 * beta * beta * sin(theta) * cos(theta) -
	0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (0.4e1 * a5 * g11 * g11 * h11
	* h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 *
	h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 *
	h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21
	* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha) * (a5 * g11 * g11
	* h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 *
	a6 * h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) *
	cos(theta) + a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) +
	a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *
	g11 * h11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 *
	g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 *
	beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta
	* pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta)
	* cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5
	* g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12
	* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22
	* beta * beta * pow(sin(theta), 0.2e1));

	H4 = -0.2e1 * beta * beta * sin(theta) * cos(theta) * (a5 * g11 * g11 * h11 *
	h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 * h12
	* h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) + a5
	* g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 * g22
	* h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11 *
	g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22 *
	beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta * beta
	* sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) +
	a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *
	g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 *
	g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 *
	beta * beta * pow(sin(theta), 0.2e1)) - (-0.2e1 * beta * pow(sin(theta), 0.2e1)
	* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (0.4e1 * a5 * g21 *
	g21 * h11 * h21 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 *
	h11 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 *
	h12 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * h12
	* g11 * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 *
	g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21
	* h21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g22 *
	h22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * h22
	* alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * beta * sin(theta) *
	cos(theta) - 0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (-0.2e1 * a5 *
	g11 * g11 * h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 *
	h11 * h21 - 0.2e1 * a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 *
	beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta
	* beta * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta *
	beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta *
	beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta
	* sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha *
	cos(theta) * sin(theta)) - (-0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha) * (0.4e1 * a5 * g11 * g11 * h11 *
	h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 *
	h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 *
	h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21
	* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta);

	J0 = -beta * cos(theta) * (-a4 * a6 * g12 * h22 - a4 * a6 * g11 * h21 + a5 * g12
	* h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * g12 * g21 * h21 * beta
	* beta * pow(cos(theta), 0.2e1) + a5 * g12 * g12 * g11 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g12 * g12 * g21 * h21 * beta * beta * sin(theta) *
	cos(theta) - a5 * g22 * g12 * g11 * h21 * beta * beta * sin(theta) * cos(theta)
	- a5 * g22 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h21 * beta
	* beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta *
	beta * sin(theta) * cos(theta) - a5 * g12 * h22 * beta * beta * sin(theta) *
	cos(theta) - a5 * g11 * h21 * beta * beta * sin(theta) * cos(theta) - a5 * g21 *
	g11 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * pow(g12, 0.3e1) *
	h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * g22 * g12 * h22 * beta *
	beta * pow(cos(theta), 0.2e1) - a5 * g21 * h21 * beta * beta * sin(theta) *
	cos(theta) + a5 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 *
	pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g11 * g11 *
	g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * g11 * g22 * h22 *
	beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * g11 * h21 * beta * beta
	* pow(cos(theta), 0.2e1) + a5 * g21 * g11 * g22 * h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) * cos(theta));

	J1 = -beta * cos(theta) * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta)
	+ 0.4e1 * a5 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *
	pow(g11, 0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 *
	g12 * h22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * g12 * g21
	* h21 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * g12 * g11 *
	h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g12 * g11 * h21
	* beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h21 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta *
	sin(theta) - 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1)
	* beta + 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1
	* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 *
	a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h21
	* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g22 * h22 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta
	* sin(theta)) + 0.2e1 * alpha * sin(theta) * (-a4 * a6 * g12 * h22 - a4 * a6 *
	g11 * h21 + a5 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 *
	g12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g12 * g12 * g11 *
	h21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g12 * g12 * g21 * h21 * beta *
	beta * sin(theta) * cos(theta) - a5 * g22 * g12 * g11 * h21 * beta * beta *
	sin(theta) * cos(theta) - a5 * g22 * h22 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11
	* g21 * h21 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h22 * beta *
	beta * sin(theta) * cos(theta) - a5 * g11 * h21 * beta * beta * sin(theta) *
	cos(theta) - a5 * g21 * g11 * g12 * h22 * beta * beta * sin(theta) * cos(theta)
	+ a5 * pow(g12, 0.3e1) * h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 *
	g22 * g12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 * h21 * beta *
	beta * sin(theta) * cos(theta) + a5 * g21 * h21 * beta * beta * pow(cos(theta),
	0.2e1) + a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5
	* g11 * g11 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * g11
	* g22 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * g11 * h21
	* beta * beta * pow(cos(theta), 0.2e1) + a5 * g21 * g11 * g22 * h22 * beta *
	beta * pow(cos(theta), 0.2e1) + a5 * g22 * h22 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) *
	cos(theta));

	J2 = -beta * cos(theta) * (- (2 * a4 * a6 * g12 * h22) -  (2 * a4 * a6 * g11 *
	h21) - 0.2e1 * a5 *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 *
	a5 * g22 *  g12 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5
	*  (g12 * g12) *  g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5
	*  (g12 * g12) * g21 *  h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5
	* g22 *  g12 *  g11 *  h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5
	* g22 *  h22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  g11 *  h21
	* beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a5 *  (g11 * g11) * g21 *  h21
	* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g12 *  h22 * beta * beta
	* sin(theta) * cos(theta) + 0.2e1 * a5 *  g11 *  h21 * beta * beta * sin(theta)
	* cos(theta) + 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 * beta * beta * sin(theta)
	* cos(theta) - 0.2e1 * a5 *   pow( g12,  3) *  h22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 *  g12 *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 *  h21 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a5 *   pow( g11,  3) *  h21 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a5 *  (g11 * g11) *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.2e1 * a5 *  (g11 * g11) * g22 *  h22 * beta * beta * sin(theta) * cos(theta) -
	0.2e1 * a5 * g21 * g21 *  g11 *  h21 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a5 * g21 *  g11 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) -
	0.2e1 * a5 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *
	(g12 * g12) * g22 *  h22 * beta * beta * sin(theta) * cos(theta) + 0.4e1 * a5 *
	(g11 * g11) *  g12 *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5
	*  (g11 * g11) * g22 *  h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 *
	a5 * g21 *  g11 *  g12 *  h22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1
	* a5 *  (g11 * g11) * g21 *  h21 * alpha * alpha * cos(theta) * sin(theta) +
	0.4e1 * a5 *   pow( g12,  3) *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a5 * g22 * g22 *  g12 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a5 *  g11 *  h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *
	g21 *  h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 *  g12 *  h22 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *   pow( g11,  3) *  h21 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *  g12 *  h22 * alpha *
	alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 *  g12 * g21 *  h21 * alpha *
	alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 *  g12 *  g11 *  h21 * alpha *
	alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *  g11 *  h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g21 *  g11 * g22 *  h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 *  h22 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 *  (g12 * g12) * g21 *  h21 * alpha * alpha * cos(theta)
	* sin(theta) + 0.4e1 * a5 * g21 *  h21 * alpha * alpha * cos(theta) * sin(theta)
	+ 0.4e1 * a5 *  (g12 * g12) *  g11 *  h21 * alpha * alpha * pow(cos(theta),
	0.2e1) + 0.4e1 * a5 * g22 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.8e1 * a5 *  (g12 * g12) * g22 *  h22 * alpha * alpha * cos(theta) * sin(theta)
	+ 0.4e1 * a5 * g21 * g21 *  g11 *  h21 * alpha * alpha * pow(sin(theta), 0.2e1))
	+ 0.2e1 * alpha * sin(theta) * (0.4e1 * a5 *  (g11 * g11) *  g12 *  h22 * alpha
	* cos(theta) * beta * sin(theta) - 0.2e1 * a5 *  (g11 * g11) * g22 *  h22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 *  (g11 * g11) * g21 *  h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 *  (g11 * g11) * g21 *  h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 *  g12 * g21 *  h21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *   pow( g12,  3) *  h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *  g12 *  h22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 *  h21 * alpha * cos(theta) *
	beta * sin(theta) + 0.4e1 * a5 *  g12 *  h22 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 *   pow( g11,  3) *  h21 * alpha * cos(theta) * beta *
	sin(theta) - 0.2e1 * a5 *  g12 *  h22 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a5 *  (g12 * g12) * g21 *  h21 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.2e1 * a5 * g22 *  g12 *  g11 *  h21 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a5 * g22 *  g12 *  g11 *  h21 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.4e1 * a5 *  g11 *  h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 *
	g11 *  h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 *  h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 *  (g12 * g12) *  g11 *  h21
	* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 *  (g12 * g12) * g21 *
	h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 *  h22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * g22 *  g12 *  h22 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g21 *  h21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 *  h21 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g21 *  g11 * g22 *  h22 * alpha * cos(theta) *
	beta * sin(theta) - 0.4e1 * a5 * g22 *  h22 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 *  (g12 * g12) * g22 *  h22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 *  (g12 * g12) * g22 *  h22 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a5 *  g11 *  h21 * beta * pow(sin(theta), 0.2e1) *
	alpha + 0.2e1 * a5 *  (g11 * g11) * g22 *  h22 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.4e1 * a5 * g21 * g21 *  g11 *  h21 * alpha * cos(theta) * beta *
	sin(theta)) + beta * cos(theta) * (- (a4 * a6 * g12 * h22) -  (a4 * a6 * g11 *
	h21) + a5 *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 *  g12
	* g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 *  (g12 * g12) *  g11 *
	h21 * beta * beta * pow(sin(theta), 0.2e1) - a5 *  (g12 * g12) * g21 *  h21 *
	beta * beta * sin(theta) * cos(theta) - a5 * g22 *  g12 *  g11 *  h21 * beta *
	beta * sin(theta) * cos(theta) - a5 * g22 *  h22 * beta * beta * sin(theta) *
	cos(theta) + a5 *  g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 *
	a5 *  (g11 * g11) * g21 *  h21 * beta * beta * sin(theta) * cos(theta) - a5 *
	g12 *  h22 * beta * beta * sin(theta) * cos(theta) - a5 *  g11 *  h21 * beta *
	beta * sin(theta) * cos(theta) - a5 * g21 *  g11 *  g12 *  h22 * beta * beta *
	sin(theta) * cos(theta) + a5 *   pow( g12,  3) *  h22 * beta * beta *
	pow(sin(theta), 0.2e1) + a5 * g22 * g22 *  g12 *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g21 *  h21 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 *   pow( g11,  3)
	*  h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 *  (g11 * g11) *  g12 *  h22
	* beta * beta * pow(sin(theta), 0.2e1) - a5 *  (g11 * g11) * g22 *  h22 * beta *
	beta * sin(theta) * cos(theta) + a5 * g21 * g21 *  g11 *  h21 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g21 *  g11 * g22 *  h22 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1)
	- 0.2e1 * a5 *  (g12 * g12) * g22 *  h22 * beta * beta * sin(theta) *
	cos(theta));

	J3 = -beta * cos(theta) * (-0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta)
	- 0.4e1 * a5 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 *
	pow(g11, 0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g12 * h22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * g12 * g21
	* h21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g22 * g12 * g11 *
	h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * g12 * g11 * h21
	* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11 * h21 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta *
	sin(theta) + 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1)
	* beta - 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1
	* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 *
	a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h21
	* beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 * g11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * h22 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta
	* sin(theta)) + 0.2e1 * alpha * sin(theta) * (-0.2e1 * a4 * a6 * g12 * h22 -
	0.2e1 * a4 * a6 * g11 * h21 - 0.2e1 * a5 * g12 * h22 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g12 * g12 * g11 * h21 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * h22 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) +
	0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) +
	0.2e1 * a5 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *
	g11 * h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g21 * g11 * g12
	* h22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * pow(g12, 0.3e1) *
	h22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * g12 * h22
	* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * h21 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta),
	0.2e1) - 0.2e1 * a5 * g11 * g11 * g12 * h22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cos(theta),
	0.2e1) - 0.2e1 * a5 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 *
	a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) + 0.4e1 * a5
	* g11 * g11 * g12 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *
	g11 * g11 * g22 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g21 * g11 * g12 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a5 *
	g11 * g11 * g21 * h21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	pow(g12, 0.3e1) * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *
	g22 * g22 * g12 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 *
	g11 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h21 *
	alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h22 * alpha * alpha
	* pow(cos(theta), 0.2e1) + 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h22 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h21 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha *
	pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h22 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g12 * g12 * g21 * h21 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g21 * h21 * alpha * alpha * cos(theta) * sin(theta) +
	0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a5 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a5 *
	g12 * g12 * g22 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g21 * g21 * g11 * h21 * alpha * alpha * pow(sin(theta), 0.2e1)) + beta *
	cos(theta) * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cos(theta) * beta *
	sin(theta) - 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * pow(cos(theta), 0.2e1)
	* beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) *
	beta + 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * pow(sin(theta), 0.2e1) *
	alpha - 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * cos(theta) * beta *
	sin(theta) + 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *
	g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * pow(g11,
	0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g12 * h22 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * g12 * g21 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * g12 * g11 * h21 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g12 * g11 * h21 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h21 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta *
	sin(theta) - 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1)
	* beta + 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1
	* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 *
	a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h21
	* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g11 * g22 * h22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g22 * h22 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta
	* sin(theta));

	K0 = -beta * sin(theta) * (a5 * pow(g21, 0.3e1) * h11 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * beta
	* beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12 * g11 * g21 * beta * beta *
	sin(theta) * cos(theta) - a5 * g11 * h11 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * h12 * beta
	* beta * pow(cos(theta), 0.2e1) + a5 * pow(g22, 0.3e1) * h12 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) +
	a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 *
	h11 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * h12 * g11 * g21 * beta
	* beta * pow(sin(theta), 0.2e1) - a5 * g12 * h12 * g21 * g21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta *
	sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 0.2e1) -
	a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12
	* g22 * g22 * beta * beta * sin(theta) * cos(theta) - a4 * a6 * h11 * g21 - a4 *
	a6 * h12 * g22);

	K1 = -beta * sin(theta) * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) * beta * sin(theta) -
	0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12
	* h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cos(theta) * beta
	* sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *
	g11 * h11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11
	* h11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 *
	h11 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g11 * h11
	* g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 *
	g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12
	* g22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha) -
	0.2e1 * alpha * cos(theta) * (a5 * pow(g21, 0.3e1) * h11 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * beta
	* beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12 * g11 * g21 * beta * beta *
	sin(theta) * cos(theta) - a5 * g11 * h11 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * h12 * beta
	* beta * pow(cos(theta), 0.2e1) + a5 * pow(g22, 0.3e1) * h12 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) +
	a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 *
	h11 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * h12 * g11 * g21 * beta
	* beta * pow(sin(theta), 0.2e1) - a5 * g12 * h12 * g21 * g21 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta *
	sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) *
	cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 0.2e1) -
	a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12
	* g22 * g22 * beta * beta * sin(theta) * cos(theta) - a4 * a6 * h11 * g21 - a4 *
	a6 * h12 * g22);

	K2 = -beta * sin(theta) * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g12 * h12 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g22
	* h12 * g11 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g11 *
	h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * beta *
	beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sin(theta) *
	cos(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.4e1
	* a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 *
	a4 * a6 * h11 * g21 - 0.2e1 * a4 * a6 * h12 * g22 + 0.4e1 * a5 * pow(g22, 0.3e1)
	* h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * pow(g21, 0.3e1) *
	h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * alpha *
	alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a5 * g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g12 * h12 * g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 *
	g11 * h11 * g21 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 *
	g22 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 *
	g12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * g21 *
	alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * g21 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * g12 * h12 * g22 * g22 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha
	* cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha *
	pow(sin(theta), 0.2e1)) - 0.2e1 * alpha * cos(theta) * (0.2e1 * a5 * g22 * h12 *
	g11 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 *
	beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 *
	alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * alpha * cos(theta) *
	beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1)
	* h11 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha
	* pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h11 * alpha * cos(theta) *
	beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cos(theta) * beta *
	sin(theta) + 0.2e1 * a5 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.2e1 * a5 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 *
	g12 * h12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12
	* h12 * g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 *
	pow(g22, 0.3e1) * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha) + beta * sin(theta) * (a5 * pow(g21,
	0.3e1) * h11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta *
	beta * sin(theta) * cos(theta) + a5 * g11 * h11 * beta * beta * pow(sin(theta),
	0.2e1) + a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12
	* g11 * g21 * beta * beta * sin(theta) * cos(theta) - a5 * g11 * h11 * beta *
	beta * sin(theta) * cos(theta) + a5 * g12 * h12 * beta * beta * pow(sin(theta),
	0.2e1) + a5 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a5 * pow(g22,
	0.3e1) * h12 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21
	* beta * beta * pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta *
	beta * sin(theta) * cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g21 * h11 * g22 * g22 * beta * beta *
	pow(cos(theta), 0.2e1) - a5 * g21 * h11 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g12 *
	h12 * g21 * g21 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11
	* g21 * g21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22
	* beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta *
	beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta *
	pow(sin(theta), 0.2e1) - a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta)
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) -
	a4 * a6 * h11 * g21 - a4 * a6 * h12 * g22);

	K3 = -beta * sin(theta) * (-0.2e1 * a5 * g22 * h12 * g11 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cos(theta) * beta *
	sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cos(theta) * beta *
	sin(theta) - 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) * beta * sin(theta) +
	0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 *
	g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12
	* h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g22 * h12 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cos(theta) * beta
	* sin(theta) + 0.2e1 * a5 * g11 * h11 * alpha * pow(cos(theta), 0.2e1) * beta -
	0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 *
	g11 * h11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11
	* h11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 *
	h11 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h11
	* g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 *
	g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g12
	* g22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g22 * h12 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha *
	cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g21 * h11 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.2e1 * a5 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha) -
	0.2e1 * alpha * cos(theta) * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta
	* pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g12 * h12 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) -
	0.2e1 * a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g22
	* h12 * g11 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g11 *
	h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * beta *
	beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * h12 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta *
	sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sin(theta) *
	cos(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) *
	cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta),
	0.2e1) + 0.2e1 * a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.4e1
	* a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 *
	a4 * a6 * h11 * g21 - 0.2e1 * a4 * a6 * h12 * g22 + 0.4e1 * a5 * pow(g22, 0.3e1)
	* h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * pow(g21, 0.3e1) *
	h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * alpha *
	alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha *
	cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha *
	pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * cos(theta) *
	sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) +
	0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) +
	0.4e1 * a5 * g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g12 * h12 * g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 *
	g11 * h11 * g21 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *
	g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 *
	g22 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 *
	g12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * g21 *
	alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * g21 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 *
	alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * g12 * h12 * g22 * g22 *
	alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha
	* cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha *
	pow(sin(theta), 0.2e1)) + beta * sin(theta) * (0.2e1 * a5 * g22 * h12 * g11 *
	g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) *
	beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta *
	sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1)
	* alpha - 0.2e1 * a5 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1
	* a5 * g22 * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 *
	h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1) * h11
	* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta *
	pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h11 * alpha * cos(theta) *
	beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cos(theta) * beta *
	sin(theta) + 0.2e1 * a5 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.2e1 * a5 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 *
	g12 * h12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12
	* h12 * g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 *
	pow(g22, 0.3e1) * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 *
	alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha);

	L0 = a4 * a6 * g11 * g21 + a4 * a6 * g12 * g22 + a4 * a6 - a5 * g22 * g22 * beta
	* beta * pow(cos(theta), 0.2e1) - a5 * g21 * g21 * beta * beta * pow(cos(theta),
	0.2e1) + a5 * g11 * g11 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - a5
	* beta * beta * pow(sin(theta), 0.2e1) - a5 * beta * beta * pow(cos(theta),
	0.2e1) - a5 * pow(g12, 0.3e1) * g22 * beta * beta * pow(sin(theta), 0.2e1) - a5
	* g12 * g12 * g11 * g21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * pow(g21,
	0.3e1) * g11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g22 * g12 * beta *
	beta * pow(cos(theta), 0.2e1) - a5 * g21 * g11 * beta * beta * pow(sin(theta),
	0.2e1) - a5 * pow(g22, 0.3e1) * g12 * beta * beta * pow(cos(theta), 0.2e1) +
	0.2e1 * a5 * g21 * g11 * beta * beta * sin(theta) * cos(theta) - a5 * g22 * g12
	* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * g11 * g12 * g22 *
	beta * beta * sin(theta) * cos(theta) - a5 * g21 * g11 * beta * beta *
	pow(cos(theta), 0.2e1) + a5 * g12 * g12 * beta * beta * sin(theta) * cos(theta)
	- a5 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * g11
	* g21 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * g12 *
	beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * beta * beta * sin(theta) *
	cos(theta) - a5 * g22 * g12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) +
	0.2e1 * a5 * g12 * g12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) + a5
	* g12 * g12 * g21 * g21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * g11
	* beta * beta * sin(theta) * cos(theta) - a5 * g11 * g11 * g12 * g22 * beta *
	beta * pow(sin(theta), 0.2e1) - a5 * pow(g11, 0.3e1) * g21 * beta * beta *
	pow(sin(theta), 0.2e1) + a5 * g21 * g21 * beta * beta * sin(theta) * cos(theta)
	+ a5 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - a5 * g21 * g11 * g22
	* g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g11 * g11 * beta * beta *
	pow(sin(theta), 0.2e1);

	L1 = 0.4e1 * a5 * pow(g22, 0.3e1) * g12 * alpha * cos(theta) * beta * sin(theta)
	+ 0.4e1 * a5 * pow(g21, 0.3e1) * g11 * alpha * cos(theta) * beta * sin(theta) +
	0.4e1 * a5 * g21 * g21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *
	g21 * g11 * g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g21
	* g11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 *
	g11 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g21 * g11 * beta *
	pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * g11 * g21 * alpha *
	cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g12, 0.3e1) * g22 * alpha *
	cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * g11 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 * g11 * beta * pow(sin(theta),
	0.2e1) * alpha + 0.4e1 * a5 * g12 * g12 * g22 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta + 0.4e1 * a5 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5
	* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g12 * g21 * g21 * alpha *
	cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g12 * alpha *
	pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g12 * beta * pow(sin(theta),
	0.2e1) * alpha - 0.2e1 * a5 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha +
	0.4e1 * a5 * g11 * g11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta +
	0.4e1 * a5 * g22 * g22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *
	g12 * g12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * g12 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * g22 *
	alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 * g11 * g22 * g22 *
	beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 * g11 * g22 * g22 *
	alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g11, 0.3e1) * g21 *
	alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * g21 * alpha *
	pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g22 * alpha * pow(cos(theta),
	0.2e1) * beta - 0.2e1 * a5 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha -
	0.4e1 * a5 * g11 * g11 * g12 * g22 * alpha * cos(theta) * beta * sin(theta) -
	0.4e1 * a5 * g11 * g11 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 *
	g12 * g12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12
	* g12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 *
	g12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha;

	L2 =  (2 * a4 * a6 * g11 * g21) +  (2 * a4 * a6 * g12 * g22) +  (2 * a4 * a6) +
	0.2e1 * a5 *  (g22 * g22) * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 *
	(g21 * g21) * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *  (g11 * g11) *
	(g22 * g22) * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * beta * beta *
	pow(sin(theta), 0.2e1) + 0.2e1 * a5 * beta * beta * pow(cos(theta), 0.2e1) +
	0.2e1 * a5 *   pow( g12,  3) *  g22 * beta * beta * pow(sin(theta), 0.2e1) +
	0.2e1 * a5 *  (g12 * g12) *  g11 *  g21 * beta * beta * pow(sin(theta), 0.2e1) +
	0.2e1 * a5 *   pow( g21,  3) *  g11 * beta * beta * pow(cos(theta), 0.2e1) +
	0.2e1 * a5 *  g22 *  g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 *
	g21 *  g11 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 *   pow( g22,  3)
	*  g12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  g21 *  g11 * beta
	* beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g22 *  g12 * beta * beta *
	pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  g21 *  g11 *  g12 *  g22 * beta * beta *
	sin(theta) * cos(theta) + 0.2e1 * a5 *  g21 *  g11 * beta * beta *
	pow(cos(theta), 0.2e1) - 0.2e1 * a5 *  (g12 * g12) * beta * beta * sin(theta) *
	cos(theta) + 0.2e1 * a5 *  (g12 * g12) * beta * beta * pow(sin(theta), 0.2e1) -
	0.4e1 * a5 *  (g11 * g11) *  (g21 * g21) * beta * beta * sin(theta) * cos(theta)
	- 0.4e1 * a5 *  g22 *  g12 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a5
	* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g22 *  g12 *  (g21 *
	g21) * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g12 * g12) *  (g22
	* g22) * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  (g12 * g12) *
	(g21 * g21) * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  (g11 * g11)
	* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  (g11 * g11) *  g12 *
	g22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 *   pow( g11,  3) *  g21
	* beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 *  (g21 * g21) * beta * beta
	* sin(theta) * cos(theta) - 0.2e1 * a5 *  (g22 * g22) * beta * beta * sin(theta)
	* cos(theta) + 0.2e1 * a5 *  g21 *  g11 *  (g22 * g22) * beta * beta *
	pow(cos(theta), 0.2e1) + 0.2e1 * a5 *  (g11 * g11) * beta * beta *
	pow(sin(theta), 0.2e1) - 0.4e1 * a5 * alpha * alpha * pow(cos(theta), 0.2e1) -
	0.4e1 * a5 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *   pow( g22,
	3) *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g21 * g21) *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a5 *  (g11 * g11) *  (g21 *
	g21) * alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g22 *  g12 *
	alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g11 * g11) * alpha *
	alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g22 * g22) * alpha * alpha *
	pow(sin(theta), 0.2e1) - 0.8e1 * a5 * alpha * alpha * cos(theta) * sin(theta) -
	0.4e1 * a5 *   pow( g21,  3) *  g11 * alpha * alpha * pow(sin(theta), 0.2e1) -
	0.4e1 * a5 *  g22 *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *
	(g12 * g12) *  g11 *  g21 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5
	*  g21 *  g11 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *   pow(
	g12,  3) *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g11 *
	g11) * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a5 *  (g12 * g12) *
	(g22 * g22) * alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g12 *
	g12) * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a5 *  g21 *  g11 *  g12
	*  g22 * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a5 *  g22 *  g12 *
	alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g21 * g21) * alpha *
	alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g21 *  g11 * alpha * alpha *
	pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g11 * g11) *  (g22 * g22) * alpha *
	alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g12 * g12) *  (g21 * g21) *
	alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *   pow( g11,  3) *  g21 *
	alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g22 * g22) * alpha *
	alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g11 * g11) *  g12 *  g22 *
	alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a5 *  g21 *  g11 * alpha *
	alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g22 *  g12 *  (g21 * g21) *
	alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g12 * g12) * alpha *
	alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g21 *  g11 *  (g22 * g22) *
	alpha * alpha * pow(sin(theta), 0.2e1);
	#endif

	A0 = a4 * a4 * g12 * g12
	+ a4 * a4 * g11 * g11
	+ a4 * a4
	+ 2.0 * a2 * pow(g11, 3.0) * g21 * beta * beta * sth * cth
	+ 2.0 * a2 * g21 * g11 * g12 * g12 * beta * beta * sth * cth
	- 2.0 * a2 * g11 * g11 * g12 * g12 * beta * beta * pow(sth, 2.0)
	- a2 * pow(g12, 4.0) * beta * beta * pow(sth, 2.0)
	- a2 * g21 * g21 * g11 * g11 * beta * beta * pow(cth, 2.0)
	+ 2.0 * a2 * g12 * g12 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * g22 * g12 * beta * beta * sth * cth
	- a2 * beta * beta * pow(cth, 2.0)
	+ 2.0 * a2 * pow(g12, 3.0) * g22 * beta * beta * sth * cth
	- a2 * pow(g11, 4.0) * beta * beta * pow(sth, 2.0)
	- 2.0 * a2 * g11 * g11 * beta * beta * pow(sth, 2.0)
	- 2.0 * a2 * g12 * g12 * beta * beta * pow(sth, 2.0)
	+ 2.0 * a2 * beta * beta * sth * cth
	- 2.0 * a2 * g21 * g11 * g22 * g12 * beta * beta * pow(cth, 2.0)
	- a2 * beta * beta * pow(sth, 2.0)
	+ 2.0 * a2 * g21 * g11 * beta * beta * sth * cth
	- a2 * g22 * g22 * g12 * g12 * beta * beta * pow(cth, 2.0)
	- 2.0 * a2 * g22 * g12 * beta * beta * pow(cth, 2.0)
	- 2.0 * a2 * g21 * g11 * beta * beta * pow(cth, 2.0)
	+ 2.0 * a2 * g22 * g12 * beta * beta * sth * cth;

	A1 = 0.4e1 * a2 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * g21 * g11 * g11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.8e1 * a2 * g21 * g11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * g22 * g12 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * g11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * g11 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g21 * g11 * beta * pow(sth, 0.2e1) * alpha
	- 0.8e1 * a2 * g11 * g11 * g12 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * pow(g11, 0.4e1) * alpha * sth * beta * cth
	- 0.8e1 * a2 * g11 * g11 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * g11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * pow(g12, 0.3e1) * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * pow(g12, 0.4e1) * alpha * sth * beta * cth
	- 0.8e1 * a2 * g12 * g12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * g11 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha;

	A2 = (2 * a4 * a4 * g12 * g12)
	+ (2 * a4 * a4 * g11 * g11)
	+ (2 * a4 * a4)
	+ 0.2e1 * a2 * pow(g12, 4) * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * pow(g11, 4) * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * (g11 * g11) * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * (g12 * g12) * beta * beta * pow(sth, 0.2e1)
	- 0.4e1 * a2 * beta * beta * sth * cth
	+ 0.2e1 * a2 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * g21 * g11 * (g12 * g12) * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g21 * (g11 * g11) * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * (g12 * g12) * beta * beta * sth * cth
	+ 0.4e1 * a2 * g21 * g11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g22 * g22 * (g12 * g12) * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * g22 * g12 * beta * beta * sth * cth
	- 0.4e1 * a2 * (g11 * g11) * beta * beta * sth * cth
	- 0.4e1 * a2 * g21 * g11 * beta * beta * sth * cth
	+ 0.4e1 * a2 * (g11 * g11) * (g12 * g12) * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * pow(g11, 3) * g21 * beta * beta * sth * cth
	- 0.4e1 * a2 * (g11 * g11) * g22 * g12 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g21 * g11 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * pow(g12, 3) * g22 * beta * beta * sth * cth
	- 0.4e1 * a2 * pow(g11, 4) * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a2 * (g11 * g11) * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a2 * pow(g12, 4) * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a2 * (g12 * g12) * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a2 * alpha * alpha * cth * sth
	- 0.4e1 * a2 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a2 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a2 * g22 * g12 * alpha * alpha * cth * sth
	- 0.4e1 * a2 * g21 * g21 * (g11 * g11) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a2 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a2 * g22 * g22 * (g12 * g12) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a2 * g21 * g11 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * (g11 * g11) * (g12 * g12) * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a2 * (g11 * g11) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g21 * g11 * g22 * g12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a2 * pow(g12, 3) * g22 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * g22 * g12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a2 * pow(g11, 3) * g21 * alpha * alpha * cth * sth
	- 0.8e1 * a2 * (g11 * g11) * g22 * g12 * alpha * alpha * cth * sth;

	B0 = -0.2e1 * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a4 * a4 * h12 * g12);

	B1 = -0.2e1 * beta * sth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 0.4e1 * alpha * cth * (a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a4 * a4 * h12 * g12);

	B2 = -0.2e1 * beta * sth * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 0.2e1 * a4 * a4 * h11 * g11
	- 0.2e1 * a4 * a4 * h12 * g12
	- 0.2e1 * a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * g11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g22 * g12 * beta * beta * pow(cth, 0.2e1))
	- 0.4e1 * alpha * cth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	+ 0.2e1 * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * pow(cth, 0.2e1)
	+ a2 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- a4 * a4 * h12 * g12);

	B3 = -0.2e1 * beta * sth * (-0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 0.4e1 * alpha * cth * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 0.2e1 * a4 * a4 * h11 * g11
	- 0.2e1 * a4 * a4 * h12 * g12
	- 0.2e1 * a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * g11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g22 * g12 * beta * beta * pow(cth, 0.2e1))
	+ 0.2e1 * beta * sth * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g22 * h12 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a2 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * h11 * alpha * sth * beta * cth
	- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.2e1 * a2 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth);

	C0 = -beta * beta * pow(sth, 0.2e1) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1));

	C1 = -beta * beta * pow(sth, 0.2e1) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- 0.4e1 * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1));

	C2 = -beta * beta * pow(sth, 0.2e1) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1))
	- 0.4e1 * beta * sth * alpha * cth * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- (-0.2e1 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(cth, 0.2e1)) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1));

	C3 = -beta * beta * pow(sth, 0.2e1) * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- 0.4e1 * beta * sth * alpha * cth * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1))
	- (-0.2e1 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(cth, 0.2e1)) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	+ 0.4e1 * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1));

	C4 = -0.2e1 * beta * beta * pow(sth, 0.2e1) * (-a4 * a4 * h12 * h12
	+ 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1))
	- 0.4e1 * beta * sth * alpha * cth * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- (-0.2e1 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(cth, 0.2e1)) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 0.2e1 * a4 * a4 * h12 * h12
	+ 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a2 * g22 * g22 * h12 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g11 * g11 * h11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a4 * a4 * h11 * h11
	+ 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cth, 0.2e1))
	+ 0.4e1 * beta * sth * alpha * cth * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth);

	E0 = 0.2e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * pow(g22, 0.3e1) * beta * beta * cth * sth
	- a3 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a3 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g11 * g21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * g22 * beta * beta * cth * sth
	- a3 * pow(g21, 0.4e1) * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * g21 * beta * beta * cth * sth
	- a3 * g12 * g12 * g22 * g22 * beta * beta * pow(sth, 0.2e1)
	- a3 * g11 * g11 * g21 * g21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ a6 * a6 * g21 * g21
	+ a6 * a6 * g22 * g22
	+ 0.2e1 * a3 * g11 * pow(g21, 0.3e1) * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ a6 * a6
	- 0.2e1 * a3 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a3 * pow(g22, 0.4e1) * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * beta * beta * cth * sth;

	E1 = -0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.8e1 * a3 * g11 * g21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * beta * pow(cth, 0.2e1) * alpha
	+ 0.8e1 * a3 * g21 * g21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * g21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * g21 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * g21 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * pow(g21, 0.3e1) * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * pow(g21, 0.3e1) * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.8e1 * a3 * g22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g21 * g21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * pow(g22, 0.4e1) * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * g21 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * beta * pow(cth, 0.2e1) * alpha;

	E2 = -0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	- 0.4e1 * a3 * g12 * g22 * beta * beta * cth * sth
	- 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * beta * beta * cth * sth
	+ 0.2e1 * a3 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a3 * g11 * g21 * beta * beta * cth * sth
	- 0.4e1 * a3 * g22 * g22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * pow(g21, 0.4e1) * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a3 * g21 * g21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * g12 * g22 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * g11 * g21 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a6 * a6 * g21 * g21
	+ 0.2e1 * a6 * a6 * g22 * g22
	- 0.4e1 * a3 * g11 * pow(g21, 0.3e1) * beta * beta * cth * sth
	- 0.4e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a6 * a6
	- 0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a3 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a3 * g11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a3 * g12 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g12 * pow(g22, 0.3e1) * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a3 * g22 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g21 * g21 * g12 * g22 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g11 * pow(g21, 0.3e1) * alpha * alpha * sth * cth
	- 0.4e1 * a3 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a3 * g21 * g21 * alpha * alpha * sth * cth
	- 0.8e1 * a3 * g11 * g21 * alpha * alpha * sth * cth
	- 0.4e1 * a3 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a3 * g11 * g21 * g22 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * pow(g22, 0.4e1) * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a3 * beta * beta * cth * sth
	- 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a3 * g21 * g21 * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a3 * alpha * alpha * sth * cth
	- 0.4e1 * a3 * pow(g22, 0.4e1) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a3 * g22 * g22 * alpha * alpha * pow(sth, 0.2e1);

	F0 = -0.2e1 * beta * cth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * pow(g22, 0.3e1) * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sth, 0.2e1)
	+ a3 * pow(g21, 0.3e1) * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sth, 0.2e1));

	F1 = -0.2e1 * beta * cth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha)
	+ 0.4e1 * alpha * sth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * pow(g22, 0.3e1) * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sth, 0.2e1)
	+ a3 * pow(g21, 0.3e1) * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sth, 0.2e1));

	F2 = -0.2e1 * beta * cth * (-(2 * a6 * a6 * h22 * g22)
	- (2 * a6 * a6 * h21 * g21)
	+ 0.2e1 * a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	- 0.2e1 * a3 * pow(g22, 3) * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g21 * h21 * (g22 * g22) * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * h21 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * (g21 * g21) * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * h21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * pow(g21, 3) * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * (g21 * g21) * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * pow(g21, 3) * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * (g21 * g21) * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * (g22 * g22) * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g11 * h21 * (g21 * g21) * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * (g22 * g22) * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * pow(g22, 3) * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g12 * h22 * (g22 * g22) * alpha * alpha * sth * cth)
	+ 0.4e1 * alpha * sth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * (g21 * g21) * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * pow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * (g22 * g22) * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * (g22 * g22) * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g21 * h21 * (g22 * g22) * alpha * sth * beta * cth
	+ 0.2e1 * a3 * (g21 * g21) * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * (g21 * g21) * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * (g21 * g21) * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * pow(g22, 3) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * (g22 * g22) * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * (g21 * g21) * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g11 * h21 * (g22 * g22) * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha)
	+ 0.2e1 * beta * cth * (-(a6 * a6 * h22 * g22)
	- (a6 * a6 * h21 * g21)
	- a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	+ a3 * pow(g22, 3) * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g21 * h21 * (g22 * g22) * beta * beta * pow(cth, 0.2e1)
	- a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * (g21 * g21) * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sth, 0.2e1)
	+ a3 * pow(g21, 3) * h21 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sth, 0.2e1));

	F3 = -0.2e1 * beta * cth * (0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * h21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g21 * h21 * beta * pow(cth, 0.2e1) * alpha
	- 0.2e1 * a3 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g11 * h21 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.2e1 * a3 * g11 * h21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha)
	+ 0.4e1 * alpha * sth * (-0.2e1 * a6 * a6 * h22 * g22
	- 0.2e1 * a6 * a6 * h21 * g21
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	- 0.2e1 * a3 * pow(g22, 0.3e1) * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g21 * h21 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g11 * h21 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g22 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g21 * g21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a3 * g21 * h21 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g12 * h22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * pow(g21, 0.3e1) * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * g12 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g22 * g22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g11 * h21 * g21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 0.8e1 * a3 * g12 * h22 * g22 * g22 * alpha * alpha * sth * cth)
	+ 0.2e1 * beta * cth * (-0.4e1 * a3 * g21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.2e1 * a3 * g11 * h21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * sth * beta * cth
	- 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.2e1 * a3 * g11 * h21 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * pow(sth, 0.2e1) * beta
	+ 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 0.2e1 * a3 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.2e1 * a3 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha);

	G0 = -beta * beta * pow(cth, 0.2e1) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G1 = -beta * beta * pow(cth, 0.2e1) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta)
	+ 0.4e1 * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G2 = -beta * beta * pow(cth, 0.2e1) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	+ 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta)
	- (-0.2e1 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(sth, 0.2e1)) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G3 = -beta * beta * pow(cth, 0.2e1) * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha)
	+ 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- (-0.2e1 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(sth, 0.2e1)) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta)
	- 0.4e1 * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G4 = -0.2e1 * beta * beta * pow(cth, 0.2e1) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth)
	+ 0.4e1 * beta * cth * alpha * sth * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha)
	- (-0.2e1 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * alpha * alpha * pow(sth, 0.2e1)) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 0.2e1 * a6 * a6 * h21 * h21
	+ 0.8e1 * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a6 * a6 * h22 * h22
	+ 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- 0.4e1 * beta * cth * alpha * sth * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cth, 0.2e1) * alpha
	+ 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sth, 0.2e1) * beta
	- 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cth, 0.2e1) * alpha
	- 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sth, 0.2e1) * beta
	+ 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sth, 0.2e1) * beta);

	H0 = -beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1));

	H1 = -beta * beta * sth * cth * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta)
	- (-0.2e1 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * alpha * pow(cth, 0.2e1) * beta) * (a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1));

	H2 = -beta * beta * sth * cth * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * alpha * pow(cth, 0.2e1) * beta) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1));

	H3 = -beta * beta * sth * cth * (0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta)
	- (-0.2e1 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * alpha * pow(cth, 0.2e1) * beta) * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta)
	- (-0.2e1 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * beta * pow(sth, 0.2e1) * alpha) * (a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1));

	H4 = -0.2e1 * beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1))
	- (-0.2e1 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * alpha * pow(cth, 0.2e1) * beta) * (0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta)
	- (-0.2e1 * beta * beta * sth * cth
	- 0.4e1 * alpha * alpha * cth * sth) * (-0.2e1 * a5 * g11 * g11 * h11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a4 * a6 * h11 * h21
	- 0.2e1 * a4 * a6 * h12 * h22
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-0.2e1 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * beta * pow(sth, 0.2e1) * alpha) * (0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta);

	J0 = -beta * cth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * pow(g12, 0.3e1) * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J1 = -beta * cth * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 0.2e1 * alpha * sth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * pow(g12, 0.3e1) * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J2 = -beta * cth * (-(2 * a4 * a6 * g12 * h22)
	- (2 * a4 * a6 * g11 * h21)
	- 0.2e1 * a5 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * (g12 * g12) * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * pow(g12, 3) * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * g22 * g12 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * pow(g11, 3) * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * (g11 * g11) * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth
	+ 0.4e1 * a5 * (g11 * g11) * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * (g11 * g11) * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * (g11 * g11) * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * pow(g12, 3) * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * pow(g11, 3) * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * (g12 * g12) * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * (g12 * g12) * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a5 * (g12 * g12) * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * alpha * pow(sth, 0.2e1))
	+ 0.2e1 * alpha * sth * (0.4e1 * a5 * (g11 * g11) * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * (g11 * g11) * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * (g11 * g11) * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * (g11 * g11) * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g12, 3) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g11, 3) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * (g12 * g12) * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * (g12 * g12) * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * (g12 * g12) * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * (g12 * g12) * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * (g12 * g12) * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * (g11 * g11) * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ beta * cth * (-(a4 * a6 * g12 * h22)
	- (a4 * a6 * g11 * h21)
	+ a5 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * (g12 * g12) * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * pow(g12, 3) * h22 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g11, 3) * h21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * (g11 * g11) * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth);

	J3 = -beta * cth * (-0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 0.2e1 * alpha * sth * (-0.2e1 * a4 * a6 * g12 * h22
	- 0.2e1 * a4 * a6 * g11 * h21
	- 0.2e1 * a5 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g12 * g12 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h21 * beta * beta * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g12 * h22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * pow(g12, 0.3e1) * h22 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * g22 * g12 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g11 * g11 * g12 * h22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g22 * h22 * beta * beta * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * g11 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.8e1 * a5 * g11 * g11 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h22 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.8e1 * a5 * g12 * g12 * g22 * h22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * alpha * pow(sth, 0.2e1))
	+ beta * cth * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * g12 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g22 * h22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * h21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth);

	K0 = -beta * sth * (a5 * pow(g21, 0.3e1) * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g22, 0.3e1) * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K1 = -beta * sth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha)
	- 0.2e1 * alpha * cth * (a5 * pow(g21, 0.3e1) * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g22, 0.3e1) * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K2 = -beta * sth * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g22 * h12 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a4 * a6 * h11 * g21
	- 0.2e1 * a4 * a6 * h12 * g22
	+ 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * pow(sth, 0.2e1))
	- 0.2e1 * alpha * cth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha)
	+ beta * sth * (a5 * pow(g21, 0.3e1) * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * pow(g22, 0.3e1) * h12 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K3 = -beta * sth * (-0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha)
	- 0.2e1 * alpha * cth * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g12 * h12 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g11 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * g22 * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * h11 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g22 * h12 * beta * beta * sth * cth
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 0.2e1 * a4 * a6 * h11 * g21
	- 0.2e1 * a4 * a6 * h12 * g22
	+ 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * h11 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.4e1 * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * alpha * pow(sth, 0.2e1)
	+ 0.4e1 * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	+ 0.8e1 * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * pow(sth, 0.2e1))
	+ beta * sth * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g12 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g12 * h12 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g12 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cth * beta * sth
	- 0.2e1 * a5 * g11 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g22 * h12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g22 * h12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * h11 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * h11 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha);

L0 = a4 * a6 * g11 * g21
	+ a4 * a6 * g12 * g22
	+ a4 * a6
	- a5 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g11 * g11 * g22 * g22 * beta * beta * sth * cth
	- a5 * beta * beta * pow(sth, 0.2e1)
	- a5 * beta * beta * pow(cth, 0.2e1)
	- a5 * pow(g12, 0.3e1) * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * g12 * g12 * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	- a5 * pow(g21, 0.3e1) * g11 * beta * beta * pow(cth, 0.2e1)
	- a5 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	- a5 * g21 * g11 * beta * beta * pow(sth, 0.2e1)
	- a5 * pow(g22, 0.3e1) * g12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * sth * cth
	- a5 * g22 * g12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * beta * beta * pow(cth, 0.2e1)
	+ a5 * g12 * g12 * beta * beta * sth * cth
	- a5 * g12 * g12 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * g11 * g11 * g21 * g21 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * sth * cth
	+ 0.2e1 * a5 * beta * beta * sth * cth
	- a5 * g22 * g12 * g21 * g21 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g12 * g12 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * g11 * beta * beta * sth * cth
	- a5 * g11 * g11 * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	- a5 * pow(g11, 0.3e1) * g21 * beta * beta * pow(sth, 0.2e1)
	+ a5 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g22 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * g22 * g22 * beta * beta * pow(cth, 0.2e1)
	- a5 * g11 * g11 * beta * beta * pow(sth, 0.2e1);

L1 = 0.4e1 * a5 * pow(g22, 0.3e1) * g12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * pow(g21, 0.3e1) * g11 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g21 * g11 * g12 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g21 * g11 * g12 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g21 * g11 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g21 * g11 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * g12 * g11 * g21 * alpha * cth * beta * sth
	- 0.4e1 * a5 * pow(g12, 0.3e1) * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g11 * g11 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g11 * g11 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g12 * g12 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g12 * g12 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g12 * g21 * g21 * alpha * cth * beta * sth
	+ 0.4e1 * a5 * g22 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.4e1 * a5 * g22 * g12 * beta * pow(sth, 0.2e1) * alpha
	- 0.2e1 * a5 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * g11 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g11 * g11 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.4e1 * a5 * g22 * g22 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g12 * g12 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * g12 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g11 * g11 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g11 * g11 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.4e1 * a5 * g21 * g11 * g22 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * pow(g11, 0.3e1) * g21 * alpha * cth * beta * sth
	+ 0.2e1 * a5 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	+ 0.2e1 * a5 * g22 * g22 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	- 0.4e1 * a5 * g11 * g11 * g12 * g22 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g11 * g11 * alpha * cth * beta * sth
	- 0.4e1 * a5 * g12 * g12 * g22 * g22 * beta * pow(sth, 0.2e1) * alpha
	+ 0.2e1 * a5 * g12 * g12 * g21 * g21 * alpha * pow(cth, 0.2e1) * beta
	- 0.2e1 * a5 * g12 * g12 * g21 * g21 * beta * pow(sth, 0.2e1) * alpha;

L2 = (2 * a4 * a6 * g11 * g21)
	+ (2 * a4 * a6 * g12 * g22)
	+ (2 * a4 * a6)
	+ 0.2e1 * a5 * (g22 * g22) * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * (g21 * g21) * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * (g11 * g11) * (g22 * g22) * beta * beta * sth * cth
	+ 0.2e1 * a5 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * pow(g12, 3) * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * (g12 * g12) * g11 * g21 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * pow(g21, 3) * g11 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * pow(g22, 3) * g12 * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a5 * g21 * g11 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * beta * beta * pow(sth, 0.2e1)
	- 0.4e1 * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * beta * beta * pow(cth, 0.2e1)
	- 0.2e1 * a5 * (g12 * g12) * beta * beta * sth * cth
	+ 0.2e1 * a5 * (g12 * g12) * beta * beta * pow(sth, 0.2e1)
	- 0.4e1 * a5 * (g11 * g11) * (g21 * g21) * beta * beta * sth * cth
	- 0.4e1 * a5 * g22 * g12 * beta * beta * sth * cth
	- 0.4e1 * a5 * beta * beta * sth * cth
	+ 0.2e1 * a5 * g22 * g12 * (g21 * g21) * beta * beta * pow(cth, 0.2e1)
	- 0.4e1 * a5 * (g12 * g12) * (g22 * g22) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g12 * g12) * (g21 * g21) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g11 * g11) * beta * beta * sth * cth
	+ 0.2e1 * a5 * (g11 * g11) * g12 * g22 * beta * beta * pow(sth, 0.2e1)
	+ 0.2e1 * a5 * pow(g11, 3) * g21 * beta * beta * pow(sth, 0.2e1)
	- 0.2e1 * a5 * (g21 * g21) * beta * beta * sth * cth
	- 0.2e1 * a5 * (g22 * g22) * beta * beta * sth * cth
	+ 0.2e1 * a5 * g21 * g11 * (g22 * g22) * beta * beta * pow(cth, 0.2e1)
	+ 0.2e1 * a5 * (g11 * g11) * beta * beta * pow(sth, 0.2e1)
	- 0.4e1 * a5 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * pow(g22, 3) * g12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * (g21 * g21) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a5 * (g11 * g11) * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g22 * g12 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * (g11 * g11) * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * (g22 * g22) * alpha * alpha * pow(sth, 0.2e1)
	- 0.8e1 * a5 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * pow(g21, 3) * g11 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * g22 * g12 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * (g12 * g12) * g11 * g21 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * g21 * g11 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * pow(g12, 3) * g22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * (g11 * g11) * alpha * alpha * cth * sth
	- 0.8e1 * a5 * (g12 * g12) * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g12 * g12) * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a5 * g21 * g11 * g12 * g22 * alpha * alpha * cth * sth
	- 0.8e1 * a5 * g22 * g12 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g21 * g11 * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * (g11 * g11) * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g12 * g12) * (g21 * g21) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * pow(g11, 3) * g21 * alpha * alpha * pow(cth, 0.2e1)
	- 0.4e1 * a5 * (g22 * g22) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * (g11 * g11) * g12 * g22 * alpha * alpha * pow(cth, 0.2e1)
	- 0.8e1 * a5 * g21 * g11 * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g22 * g12 * (g21 * g21) * alpha * alpha * pow(sth, 0.2e1)
	- 0.4e1 * a5 * (g12 * g12) * alpha * alpha * cth * sth
	- 0.4e1 * a5 * g21 * g11 * (g22 * g22) * alpha * alpha * pow(sth, 0.2e1);

	//% t-dependent stuff:

	//%pose_from_point_tangents_2_fn_t(t);

	//% Not used - just to check magnitute of poly.
	//double int_coef[] = {A0, A1, A2, B0, B1, B2, B3, C0, C1, C2, C3, C4, E0, E1, E2, F0, F1, F2, F3, G0, G1, G2, G3, G4, H0, H1, H2, H3, H4, J0, J1, J2, J3, K0, K1, K2, K3, L0, L1, L2};
	//common::norm(int_coef);

}

}


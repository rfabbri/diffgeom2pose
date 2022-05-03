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

	T V[3], buf[3];
	T a1, a2, a3, a4, a5, a6;

	const T g11 = gama1[0], g12 = gama1[1],
	        g21 = gama2[0], g22 = gama2[1];

  const T g11_3 = g11*g11*g11;
  const T g11_4 = g11_3*g11;
  
  const T g12_3 = g12*g12*g12;
  const T g12_4 = g12_3*g12;

	const T h11 = tgt1[0],  h12 = tgt1[1],
	        h21 = tgt2[0],  h22 = tgt2[1];

	vec1vec2_3el_sub(Gama1, Gama2, V);
	vec_3el_wise_mult2(V, V, buf);       a1 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt1, Tgt1, buf); a2 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt2, Tgt2, buf); a3 = vec_3el_sum(buf);
	vec_3el_wise_mult2(V, Tgt1, buf);    a4 = vec_3el_sum(buf);
	vec_3el_wise_mult2(Tgt1, Tgt2, buf); a5 = vec_3el_sum(buf);
	vec_3el_wise_mult2(V, Tgt2, buf);    a6 = vec_3el_sum(buf);

	const T t4 = g11 * g11,
	        t6 = g21 * g21,
	        t7 = g22 * g22;
	const T t11 = 2 * (1 + g11 * g21 + g12 * g22) / (t4 + g12 * g12  - t6 - t7);

	theta = 0.5*atan(t11);
	if (theta < 0) theta += PI / 2;
	sth = sin(theta);
  const T sth2 = sth*sth;
  const T sth3 = sth2*sth;
	cth = cos(theta);
  const T cth2 = cth*cth;
  const T cth3 = cth2*cth;
  

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

  const T
	t1 = sin(2*theta),
	t2 = cos(2*theta),
	t5 = g22 * g22,
	t8 = g21 * g21,
	t14 = g12 * g12,
	//double t15 = t1 * t1; // double-checked: not used in matlab
	t21 = g11 * g11;

	const T den1 = 2*t1*(g11*g21 +g12*g22 + 1) + t2*(t21 + t14 - t8 - t5),
	        den2 = t21 + t14 + t8 + t5 + 2;

  const T t25 = -2*a1 / (den1 - den2);
	beta = sqrt(t25);

	const T t24 = 2*a1 / (den1 + den2);
	alpha = sqrt(t24);

	//% Coefficient code adapted from Maple ::: can be further cleaned up but works

	A0 = a4 * a4 * g12 * g12
	+ a4 * a4 * g11 * g11
	+ a4 * a4
	+ 2.0 * a2 * g11_3 * g21 * beta * beta * sth * cth
	+ 2.0 * a2 * g21 * g11 * g12 * g12 * beta * beta * sth * cth
	- 2.0 * a2 * g11 * g11 * g12 * g12 * beta * beta * sth2
	- a2 * g12_4 * beta * beta * sth2
	- a2 * g21 * g21 * g11 * g11 * beta * beta * cth2
	+ 2.0 * a2 * g12 * g12 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * beta * beta * sth * cth
	+ 2.0 * a2 * g11 * g11 * g22 * g12 * beta * beta * sth * cth
	- a2 * beta * beta * cth2
	+ 2.0 * a2 * g12_3 * g22 * beta * beta * sth * cth
	- a2 * g11_4 * beta * beta * sth2
	- 2.0 * a2 * g11 * g11 * beta * beta * sth2
	- 2.0 * a2 * g12 * g12 * beta * beta * sth2
	+ 2.0 * a2 * beta * beta * sth * cth
	- 2.0 * a2 * g21 * g11 * g22 * g12 * beta * beta * cth2
	- a2 * beta * beta * sth2
	+ 2.0 * a2 * g21 * g11 * beta * beta * sth * cth
	- a2 * g22 * g22 * g12 * g12 * beta * beta * cth2
	- 2.0 * a2 * g22 * g12 * beta * beta * cth2
	- 2.0 * a2 * g21 * g11 * beta * beta * cth2
	+ 2.0 * a2 * g22 * g12 * beta * beta * sth * cth;

	A1 = 4. * a2 * alpha * cth2 * beta
	- 4. * a2 * beta * sth2 * alpha
	+ 4. * a2 * g21 * g21 * g11 * g11 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * g12 * alpha * cth2 * beta
	+ 8. * a2 * g21 * g11 * alpha * sth * beta * cth
	- 4. * a2 * g12 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * g12 * alpha * cth2 * beta
	- 4. * a2 * g22 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * g22 * g12 * g12 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g11 * g12 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g21 * g11 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * beta * sth2 * alpha
	+ 4. * a2 * g21 * g11 * g12 * g12 * alpha * cth2 * beta
	- 4. * a2 * g21 * g11 * beta * sth2 * alpha
	- 8. * a2 * g11 * g11 * g12 * g12 * alpha * sth * beta * cth
	- 4. * a2 * g11_4 * alpha * sth * beta * cth
	- 8. * a2 * g11 * g11 * alpha * sth * beta * cth
	+ 8. * a2 * g21 * g11 * g22 * g12 * alpha * sth * beta * cth
	+ 4. * a2 * g12_3 * g22 * alpha * cth2 * beta
	- 4. * a2 * g12_3 * g22 * beta * sth2 * alpha
	- 4. * a2 * g12_4 * alpha * sth * beta * cth
	- 8. * a2 * g12 * g12 * alpha * sth * beta * cth
	+ 8. * a2 * g22 * g12 * alpha * sth * beta * cth
	+ 4. * a2 * g11_3 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11_3 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g11 * g11 * g22 * g12 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * g22 * g12 * beta * sth2 * alpha;

	A2 = (2 * a4 * a4 * g12 * g12)
	+ (2 * a4 * a4 * g11 * g11)
	+ (2 * a4 * a4)
	+ 2. * a2 * g12_4 * beta * beta * sth2
	+ 2. * a2 * g11_4 * beta * beta * sth2
	+ 4. * a2 * (g11 * g11) * beta * beta * sth2
	+ 4. * a2 * (g12 * g12) * beta * beta * sth2
	- 4. * a2 * beta * beta * sth * cth
	+ 2. * a2 * beta * beta * sth2
	+ 2. * a2 * beta * beta * cth2
	- 4. * a2 * g21 * g11 * (g12 * g12) * beta * beta * sth * cth
	+ 2. * a2 * g21 * g21 * (g11 * g11) * beta * beta * cth2
	- 4. * a2 * (g12 * g12) * beta * beta * sth * cth
	+ 4. * a2 * g21 * g11 * beta * beta * cth2
	+ 2. * a2 * g22 * g22 * (g12 * g12) * beta * beta * cth2
	- 4. * a2 * g22 * g12 * beta * beta * sth * cth
	- 4. * a2 * (g11 * g11) * beta * beta * sth * cth
	- 4. * a2 * g21 * g11 * beta * beta * sth * cth
	+ 4. * a2 * (g11 * g11) * (g12 * g12) * beta * beta * sth2
	+ 4. * a2 * g22 * g12 * beta * beta * cth2
	- 4. * a2 * g11_3 * g21 * beta * beta * sth * cth
	- 4. * a2 * (g11 * g11) * g22 * g12 * beta * beta * sth * cth
	+ 4. * a2 * g21 * g11 * g22 * g12 * beta * beta * cth2
	- 4. * a2 * g12_3 * g22 * beta * beta * sth * cth
	- 4. * a2 * g11_4 * alpha * alpha * cth2
	- 8. * a2 * (g11 * g11) * alpha * alpha * cth2
	- 4. * a2 * g12_4 * alpha * alpha * cth2
	- 8. * a2 * (g12 * g12) * alpha * alpha * cth2
	- 8. * a2 * alpha * alpha * cth * sth
	- 4. * a2 * alpha * alpha * cth2
	- 4. * a2 * alpha * alpha * sth2
	- 8. * a2 * g22 * g12 * alpha * alpha * cth * sth
	- 4. * a2 * g21 * g21 * (g11 * g11) * alpha * alpha * sth2
	- 8. * a2 * (g12 * g12) * alpha * alpha * cth * sth
	- 8. * a2 * g21 * g11 * alpha * alpha * sth2
	- 4. * a2 * g22 * g22 * (g12 * g12) * alpha * alpha * sth2
	- 8. * a2 * g21 * g11 * alpha * alpha * cth * sth
	- 8. * a2 * (g11 * g11) * (g12 * g12) * alpha * alpha * cth2
	- 8. * a2 * (g11 * g11) * alpha * alpha * cth * sth
	- 8. * a2 * g21 * g11 * (g12 * g12) * alpha * alpha * cth * sth
	- 8. * a2 * g21 * g11 * g22 * g12 * alpha * alpha * sth2
	- 8. * a2 * g12_3 * g22 * alpha * alpha * cth * sth
	- 8. * a2 * g22 * g12 * alpha * alpha * sth2
	- 8. * a2 * g11_3 * g21 * alpha * alpha * cth * sth
	- 8. * a2 * (g11 * g11) * g22 * g12 * alpha * alpha * cth * sth;

	B0 = -2. * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * cth2
	+ a2 * g12_3 * h12 * beta * beta * sth2
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * cth2
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * sth2
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * cth2
	+ a2 * g22 * h12 * beta * beta * cth2
	+ a2 * g12 * h12 * beta * beta * sth2
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 2. * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11_3 * h11 * beta * beta * sth2
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * cth2
	+ a2 * g21 * h11 * beta * beta * cth2
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * sth2
	- a4 * a4 * h12 * g12);

	B1 = -2. * beta * sth * (2. * a2 * g11 * g11 * g22 * h12 * beta * sth2 * alpha
	- 2. * a2 * g11 * g11 * g22 * h12 * alpha * cth2 * beta
	- 2. * a2 * g21 * h11 * g12 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g21 * g11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g12 * g12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g12 * g12 * h12 * g22 * beta * sth2 * alpha
	+ 2. * a2 * g21 * h11 * g12 * g12 * beta * sth2 * alpha
	- 2. * a2 * g11 * h11 * g22 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g11 * h11 * g22 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g22 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g22 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g12 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * h12 * alpha * sth * beta * cth
	- 2. * a2 * g11 * h11 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * h11 * g21 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * g21 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g12_3 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11_3 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * h11 * alpha * cth2 * beta
	+ 2. * a2 * g21 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * g11 * g12 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 2. * a2 * g11 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 4. * alpha * cth * (a2 * g21 * g11 * g22 * h12 * beta * beta * cth2
	+ a2 * g12_3 * h12 * beta * beta * sth2
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * cth2
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * sth2
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * cth2
	+ a2 * g22 * h12 * beta * beta * cth2
	+ a2 * g12 * h12 * beta * beta * sth2
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 2. * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11_3 * h11 * beta * beta * sth2
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * cth2
	+ a2 * g21 * h11 * beta * beta * cth2
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * sth2
	- a4 * a4 * h12 * g12);

	B2 = -2. * beta * sth * (4. * a2 * g21 * g21 * h11 * g11 * alpha * alpha * sth2
	+ 4. * a2 * g11 * g11 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g11 * h11 * alpha * alpha * cth2
	+ 4. * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 8. * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 4. * a2 * g12_3 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g11_3 * h11 * alpha * alpha * cth2
	+ 4. * a2 * g21 * h11 * g22 * g12 * alpha * alpha * sth2
	+ 4. * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g12 * g12 * alpha * alpha * cth2
	+ 4. * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 4. * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g22 * g22 * h12 * g12 * alpha * alpha * sth2
	+ 4. * a2 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a2 * g21 * g11 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g21 * h11 * alpha * alpha * sth2
	+ 4. * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 8. * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 2. * a4 * a4 * h11 * g11
	- 2. * a4 * a4 * h12 * g12
	- 2. * a2 * g22 * g22 * h12 * g12 * beta * beta * cth2
	+ 2. * a2 * g12 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g22 * h12 * beta * beta * sth * cth
	- 2. * a2 * g21 * g21 * h11 * g11 * beta * beta * cth2
	+ 4. * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 2. * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 2. * a2 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * g12 * g12 * beta * beta * sth2
	- 2. * a2 * g11 * g11 * g12 * h12 * beta * beta * sth2
	+ 2. * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 2. * a2 * g12 * h12 * beta * beta * sth2
	+ 2. * a2 * g11 * h11 * beta * beta * sth * cth
	- 2. * a2 * g11 * h11 * beta * beta * sth2
	- 2. * a2 * g11_3 * h11 * beta * beta * sth2
	+ 4. * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 2. * a2 * g21 * g11 * g22 * h12 * beta * beta * cth2
	+ 2. * a2 * g21 * h11 * beta * beta * sth * cth
	- 2. * a2 * g12_3 * h12 * beta * beta * sth2
	- 2. * a2 * g21 * h11 * g22 * g12 * beta * beta * cth2)
	- 4. * alpha * cth * (2. * a2 * g11 * g11 * g22 * h12 * beta * sth2 * alpha
	- 2. * a2 * g11 * g11 * g22 * h12 * alpha * cth2 * beta
	- 2. * a2 * g21 * h11 * g12 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g21 * g11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g12 * g12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g12 * g12 * h12 * g22 * beta * sth2 * alpha
	+ 2. * a2 * g21 * h11 * g12 * g12 * beta * sth2 * alpha
	- 2. * a2 * g11 * h11 * g22 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g11 * h11 * g22 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g22 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g22 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g12 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * h12 * alpha * sth * beta * cth
	- 2. * a2 * g11 * h11 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * h11 * g21 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * g21 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g12_3 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11_3 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * h11 * alpha * cth2 * beta
	+ 2. * a2 * g21 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * g11 * g12 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 2. * a2 * g11 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	+ 2. * beta * sth * (a2 * g21 * g11 * g22 * h12 * beta * beta * cth2
	+ a2 * g12_3 * h12 * beta * beta * sth2
	+ a2 * g21 * h11 * g22 * g12 * beta * beta * cth2
	+ a2 * g11 * h11 * g12 * g12 * beta * beta * sth2
	- a2 * g11 * h11 * beta * beta * sth * cth
	- a4 * a4 * h11 * g11
	+ a2 * g11 * h11 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * g12 * beta * beta * cth2
	+ a2 * g22 * h12 * beta * beta * cth2
	+ a2 * g12 * h12 * beta * beta * sth2
	- a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- a2 * g12 * h12 * beta * beta * sth * cth
	- a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- a2 * g21 * h11 * beta * beta * sth * cth
	- a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	- 2. * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- a2 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11_3 * h11 * beta * beta * sth2
	+ a2 * g21 * g21 * h11 * g11 * beta * beta * cth2
	+ a2 * g21 * h11 * beta * beta * cth2
	- a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * g12 * h12 * beta * beta * sth2
	- a4 * a4 * h12 * g12);

	B3 = -2. * beta * sth * (-2. * a2 * g11 * g11 * g22 * h12 * beta * sth2 * alpha
	+ 2. * a2 * g11 * g11 * g22 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g21 * h11 * g12 * g12 * alpha * cth2 * beta
	- 2. * a2 * g21 * g11 * g12 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g12 * g12 * h12 * g22 * alpha * cth2 * beta
	- 4. * a2 * g12 * g12 * h12 * g22 * beta * sth2 * alpha
	- 2. * a2 * g21 * h11 * g12 * g12 * beta * sth2 * alpha
	+ 2. * a2 * g11 * h11 * g22 * g12 * alpha * cth2 * beta
	- 2. * a2 * g11 * h11 * g22 * g12 * beta * sth2 * alpha
	- 4. * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	+ 2. * a2 * g22 * h12 * alpha * cth2 * beta
	- 2. * a2 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	+ 2. * a2 * g12 * h12 * alpha * cth2 * beta
	- 2. * a2 * g12 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * h12 * alpha * sth * beta * cth
	+ 2. * a2 * g11 * h11 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	- 4. * a2 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g12_3 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g11_3 * h11 * alpha * sth * beta * cth
	+ 2. * a2 * g21 * h11 * alpha * cth2 * beta
	- 2. * a2 * g21 * h11 * beta * sth2 * alpha
	+ 4. * a2 * g21 * h11 * alpha * sth * beta * cth
	+ 2. * a2 * g21 * g11 * g12 * h12 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	- 2. * a2 * g11 * h11 * beta * sth2 * alpha
	+ 4. * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth)
	- 4. * alpha * cth * (4. * a2 * g21 * g21 * h11 * g11 * alpha * alpha * sth2
	+ 4. * a2 * g11 * g11 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g11 * h11 * alpha * alpha * cth2
	+ 4. * a2 * g11 * g11 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * alpha * alpha * cth * sth
	+ 8. * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cth * sth
	+ 4. * a2 * g12_3 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a2 * g11_3 * h11 * alpha * alpha * cth2
	+ 4. * a2 * g21 * h11 * g22 * g12 * alpha * alpha * sth2
	+ 4. * a2 * g21 * h11 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g12 * g12 * alpha * alpha * cth2
	+ 4. * a2 * g21 * h11 * g12 * g12 * alpha * alpha * cth * sth
	+ 4. * a2 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g22 * g22 * h12 * g12 * alpha * alpha * sth2
	+ 4. * a2 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a2 * g21 * g11 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a2 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g21 * h11 * alpha * alpha * sth2
	+ 4. * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cth * sth
	+ 8. * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cth * sth
	- 2. * a4 * a4 * h11 * g11
	- 2. * a4 * a4 * h12 * g12
	- 2. * a2 * g22 * g22 * h12 * g12 * beta * beta * cth2
	+ 2. * a2 * g12 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g22 * h12 * beta * beta * sth * cth
	- 2. * a2 * g21 * g21 * h11 * g11 * beta * beta * cth2
	+ 4. * a2 * g12 * g12 * h12 * g22 * beta * beta * sth * cth
	+ 2. * a2 * g21 * g11 * g12 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g21 * h11 * g12 * g12 * beta * beta * sth * cth
	- 2. * a2 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * g12 * g12 * beta * beta * sth2
	- 2. * a2 * g11 * g11 * g12 * h12 * beta * beta * sth2
	+ 2. * a2 * g11 * g11 * g22 * h12 * beta * beta * sth * cth
	+ 2. * a2 * g11 * h11 * g22 * g12 * beta * beta * sth * cth
	- 2. * a2 * g12 * h12 * beta * beta * sth2
	+ 2. * a2 * g11 * h11 * beta * beta * sth * cth
	- 2. * a2 * g11 * h11 * beta * beta * sth2
	- 2. * a2 * g11_3 * h11 * beta * beta * sth2
	+ 4. * a2 * g11 * g11 * h11 * g21 * beta * beta * sth * cth
	- 2. * a2 * g21 * g11 * g22 * h12 * beta * beta * cth2
	+ 2. * a2 * g21 * h11 * beta * beta * sth * cth
	- 2. * a2 * g12_3 * h12 * beta * beta * sth2
	- 2. * a2 * g21 * h11 * g22 * g12 * beta * beta * cth2)
	+ 2. * beta * sth * (2. * a2 * g11 * g11 * g22 * h12 * beta * sth2 * alpha
	- 2. * a2 * g11 * g11 * g22 * h12 * alpha * cth2 * beta
	- 2. * a2 * g21 * h11 * g12 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g21 * g11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g12 * g12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g12 * g12 * h12 * g22 * beta * sth2 * alpha
	+ 2. * a2 * g21 * h11 * g12 * g12 * beta * sth2 * alpha
	- 2. * a2 * g11 * h11 * g22 * g12 * alpha * cth2 * beta
	+ 2. * a2 * g11 * h11 * g22 * g12 * beta * sth2 * alpha
	+ 4. * a2 * g11 * h11 * g12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g22 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g22 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * g22 * h12 * g12 * alpha * sth * beta * cth
	- 2. * a2 * g12 * h12 * alpha * cth2 * beta
	+ 2. * a2 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g22 * h12 * alpha * sth * beta * cth
	- 2. * a2 * g11 * h11 * alpha * cth2 * beta
	- 4. * a2 * g11 * g11 * h11 * g21 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * g21 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * g22 * g12 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g12_3 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11_3 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * h11 * alpha * cth2 * beta
	+ 2. * a2 * g21 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * h11 * alpha * sth * beta * cth
	- 2. * a2 * g21 * g11 * g12 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * g12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g11 * g22 * h12 * alpha * sth * beta * cth
	+ 2. * a2 * g11 * h11 * beta * sth2 * alpha
	- 4. * a2 * g21 * g21 * h11 * g11 * alpha * sth * beta * cth);

	C0 = -beta * beta * sth2 * (-a4 * a4 * h12 * h12
	+ 2. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 2. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	+ 2. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	- 2. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * sth2);

	C1 = -beta * beta * sth2 * (8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- 4. * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 2. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 2. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	+ 2. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	- 2. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * sth2);

	C2 = -beta * beta * sth2 * (-4. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 2. * a4 * a4 * h12 * h12
	+ 4. * a2 * g21 * g21 * h11 * h11 * alpha * alpha * sth2
	+ 8. * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * alpha * cth2
	- 2. * a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	- 4. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	+ 8. * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 8. * a2 * g11 * h11 * g12 * h12 * alpha * alpha * cth2
	+ 8. * a2 * g21 * h11 * g22 * h12 * alpha * alpha * sth2
	- 2. * a2 * g11 * g11 * h11 * h11 * beta * beta * sth2
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * alpha * cth2
	- 2. * a4 * a4 * h11 * h11
	+ 4. * a2 * g22 * g22 * h12 * h12 * alpha * alpha * sth2
	+ 8. * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 8. * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	- 2. * a2 * g21 * g21 * h11 * h11 * beta * beta * cth2)
	- 4. * beta * sth * alpha * cth * (8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	- (-2. * beta * beta * sth2
	+ 4. * alpha * alpha * cth2) * (-a4 * a4 * h12 * h12
	+ 2. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 2. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	+ 2. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	- 2. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * sth2);

	C3 = -beta * beta * sth2 * (4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	- 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- 4. * beta * sth * alpha * cth * (-4. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 2. * a4 * a4 * h12 * h12
	+ 4. * a2 * g21 * g21 * h11 * h11 * alpha * alpha * sth2
	+ 8. * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * alpha * cth2
	- 2. * a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	- 4. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	+ 8. * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 8. * a2 * g11 * h11 * g12 * h12 * alpha * alpha * cth2
	+ 8. * a2 * g21 * h11 * g22 * h12 * alpha * alpha * sth2
	- 2. * a2 * g11 * g11 * h11 * h11 * beta * beta * sth2
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * alpha * cth2
	- 2. * a4 * a4 * h11 * h11
	+ 4. * a2 * g22 * g22 * h12 * h12 * alpha * alpha * sth2
	+ 8. * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 8. * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	- 2. * a2 * g21 * g21 * h11 * h11 * beta * beta * cth2)
	- (-2. * beta * beta * sth2
	+ 4. * alpha * alpha * cth2) * (8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth)
	+ 4. * beta * sth * alpha * cth * (-a4 * a4 * h12 * h12
	+ 2. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 2. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	+ 2. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	- 2. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * sth2);

	C4 = -2. * beta * beta * sth2 * (-a4 * a4 * h12 * h12
	+ 2. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	- 2. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- a4 * a4 * h11 * h11
	- 2. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ a2 * g21 * g21 * h11 * h11 * beta * beta * cth2
	- 2. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	+ a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	+ a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	+ 2. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	- 2. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ a2 * g11 * g11 * h11 * h11 * beta * beta * sth2)
	- 4. * beta * sth * alpha * cth * (4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	+ 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	- 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth)
	- (-2. * beta * beta * sth2
	+ 4. * alpha * alpha * cth2) * (-4. * a2 * g11 * h11 * g12 * h12 * beta * beta * sth2
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * beta * sth * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * beta * sth * cth
	- 2. * a4 * a4 * h12 * h12
	+ 4. * a2 * g21 * g21 * h11 * h11 * alpha * alpha * sth2
	+ 8. * a2 * g21 * h11 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * alpha * cth2
	- 2. * a2 * g22 * g22 * h12 * h12 * beta * beta * cth2
	- 4. * a2 * g21 * h11 * g22 * h12 * beta * beta * cth2
	+ 8. * a2 * g11 * h11 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * beta * sth * cth
	+ 8. * a2 * g11 * h11 * g12 * h12 * alpha * alpha * cth2
	+ 8. * a2 * g21 * h11 * g22 * h12 * alpha * alpha * sth2
	- 2. * a2 * g11 * g11 * h11 * h11 * beta * beta * sth2
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * alpha * cth2
	- 2. * a4 * a4 * h11 * h11
	+ 4. * a2 * g22 * g22 * h12 * h12 * alpha * alpha * sth2
	+ 8. * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cth * sth
	+ 8. * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cth * sth
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * beta * sth * cth
	- 2. * a2 * g12 * g12 * h12 * h12 * beta * beta * sth2
	- 2. * a2 * g21 * g21 * h11 * h11 * beta * beta * cth2)
	+ 4. * beta * sth * alpha * cth * (8. * a2 * g11 * h11 * g12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g22 * g22 * h12 * h12 * alpha * sth * beta * cth
	- 4. * a2 * g21 * g21 * h11 * h11 * alpha * sth * beta * cth
	- 4. * a2 * g21 * h11 * g12 * h12 * alpha * cth2 * beta
	- 4. * a2 * g12 * h12 * h12 * g22 * alpha * cth2 * beta
	+ 4. * a2 * g21 * h11 * g12 * h12 * beta * sth2 * alpha
	- 4. * a2 * g11 * h11 * h11 * g21 * alpha * cth2 * beta
	- 4. * a2 * g11 * h11 * g22 * h12 * alpha * cth2 * beta
	+ 4. * a2 * g11 * g11 * h11 * h11 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * g22 * h12 * beta * sth2 * alpha
	+ 4. * a2 * g12 * h12 * h12 * g22 * beta * sth2 * alpha
	- 8. * a2 * g21 * h11 * g22 * h12 * alpha * sth * beta * cth
	+ 4. * a2 * g11 * h11 * h11 * g21 * beta * sth2 * alpha
	+ 4. * a2 * g12 * g12 * h12 * h12 * alpha * sth * beta * cth);

	E0 = 2. * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	+ 2. * a3 * g12 * g22 * beta * beta * cth * sth
	+ 2. * a3 * g12 * intpow(g22, 3) * beta * beta * cth * sth
	- a3 * beta * beta * sth2
	- 2. * a3 * g12 * g22 * beta * beta * sth2
	- 2. * a3 * g11 * g21 * beta * beta * sth2
	- a3 * beta * beta * cth2
	+ 2. * a3 * g11 * g21 * beta * beta * cth * sth
	+ 2. * a3 * g22 * g22 * beta * beta * cth * sth
	- a3 * intpow(g21, 4) * beta * beta * cth2
	+ 2. * a3 * g21 * g21 * beta * beta * cth * sth
	- a3 * g12 * g12 * g22 * g22 * beta * beta * sth2
	- a3 * g11 * g11 * g21 * g21 * beta * beta * sth2
	- 2. * a3 * g21 * g21 * g22 * g22 * beta * beta * cth2
	+ a6 * a6 * g21 * g21
	+ a6 * a6 * g22 * g22
	+ 2. * a3 * g11 * intpow(g21, 3) * beta * beta * cth * sth
	+ 2. * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g11 * g21 * g12 * g22 * beta * beta * sth2
	+ a6 * a6
	- 2. * a3 * g21 * g21 * beta * beta * cth2
	- 2. * a3 * g22 * g22 * beta * beta * cth2
	- a3 * intpow(g22, 4) * beta * beta * cth2
	+ 2. * a3 * beta * beta * cth * sth;

	E1 = -4. * a3 * g11 * g11 * g21 * g21 * alpha * sth * beta * cth
	+ 8. * a3 * g21 * g21 * g22 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g12 * g12 * g22 * g22 * alpha * sth * beta * cth
	+ 4. * a3 * g22 * g22 * beta * cth2 * alpha
	- 4. * a3 * g22 * g22 * alpha * sth2 * beta
	- 8. * a3 * g11 * g21 * alpha * sth * beta * cth
	- 4. * a3 * g12 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g12 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g21 * g21 * g12 * g22 * beta * cth2 * alpha
	- 4. * a3 * g12 * intpow(g22, 3) * alpha * sth2 * beta
	+ 4. * a3 * g12 * intpow(g22, 3) * beta * cth2 * alpha
	+ 8. * a3 * g21 * g21 * alpha * sth * beta * cth
	+ 4. * a3 * intpow(g21, 4) * alpha * sth * beta * cth
	- 4. * a3 * g11 * g21 * alpha * sth2 * beta
	+ 4. * a3 * g11 * g21 * beta * cth2 * alpha
	- 8. * a3 * g12 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * g12 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * g21 * g22 * g22 * beta * cth2 * alpha
	- 8. * a3 * g11 * g21 * g12 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * intpow(g21, 3) * alpha * sth2 * beta
	+ 4. * a3 * g11 * intpow(g21, 3) * beta * cth2 * alpha
	+ 4. * a3 * g21 * g21 * beta * cth2 * alpha
	+ 8. * a3 * g22 * g22 * alpha * sth * beta * cth
	- 4. * a3 * alpha * sth2 * beta
	- 4. * a3 * g21 * g21 * alpha * sth2 * beta
	+ 4. * a3 * intpow(g22, 4) * alpha * sth * beta * cth
	- 4. * a3 * g11 * g21 * g22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * beta * cth2 * alpha;

	E2 = -4. * a3 * g21 * g21 * g12 * g22 * beta * beta * cth * sth
	- 4. * a3 * g12 * g22 * beta * beta * cth * sth
	- 4. * a3 * g12 * intpow(g22, 3) * beta * beta * cth * sth
	+ 2. * a3 * beta * beta * sth2
	+ 4. * a3 * g12 * g22 * beta * beta * sth2
	+ 4. * a3 * g11 * g21 * beta * beta * sth2
	+ 2. * a3 * beta * beta * cth2
	- 4. * a3 * g11 * g21 * beta * beta * cth * sth
	- 4. * a3 * g22 * g22 * beta * beta * cth * sth
	+ 2. * a3 * intpow(g21, 4) * beta * beta * cth2
	- 4. * a3 * g21 * g21 * beta * beta * cth * sth
	+ 2. * a3 * g12 * g12 * g22 * g22 * beta * beta * sth2
	+ 2. * a3 * g11 * g11 * g21 * g21 * beta * beta * sth2
	+ 4. * a3 * g21 * g21 * g22 * g22 * beta * beta * cth2
	+ 2. * a6 * a6 * g21 * g21
	+ 2. * a6 * a6 * g22 * g22
	- 4. * a3 * g11 * intpow(g21, 3) * beta * beta * cth * sth
	- 4. * a3 * g11 * g21 * g22 * g22 * beta * beta * cth * sth
	+ 4. * a3 * g11 * g21 * g12 * g22 * beta * beta * sth2
	+ 2. * a6 * a6
	- 4. * a3 * g11 * g11 * g21 * g21 * alpha * alpha * cth2
	- 8. * a3 * g12 * g22 * alpha * alpha * cth2
	- 8. * a3 * g11 * g21 * alpha * alpha * cth2
	- 8. * a3 * g12 * g22 * alpha * alpha * sth * cth
	- 8. * a3 * g12 * intpow(g22, 3) * alpha * alpha * sth * cth
	- 8. * a3 * g21 * g21 * g22 * g22 * alpha * alpha * sth2
	- 8. * a3 * g22 * g22 * alpha * alpha * sth * cth
	- 8. * a3 * g21 * g21 * g12 * g22 * alpha * alpha * sth * cth
	- 8. * a3 * g11 * intpow(g21, 3) * alpha * alpha * sth * cth
	- 4. * a3 * alpha * alpha * sth2
	- 8. * a3 * g11 * g21 * g12 * g22 * alpha * alpha * cth2
	- 4. * a3 * g12 * g12 * g22 * g22 * alpha * alpha * cth2
	- 8. * a3 * g21 * g21 * alpha * alpha * sth * cth
	- 8. * a3 * g11 * g21 * alpha * alpha * sth * cth
	- 4. * a3 * alpha * alpha * cth2
	- 8. * a3 * g11 * g21 * g22 * g22 * alpha * alpha * sth * cth
	+ 4. * a3 * g21 * g21 * beta * beta * cth2
	+ 4. * a3 * g22 * g22 * beta * beta * cth2
	+ 2. * a3 * intpow(g22, 4) * beta * beta * cth2
	- 4. * a3 * beta * beta * cth * sth
	- 4. * a3 * intpow(g21, 4) * alpha * alpha * sth2
	- 8. * a3 * g21 * g21 * alpha * alpha * sth2
	- 8. * a3 * alpha * alpha * sth * cth
	- 4. * a3 * intpow(g22, 4) * alpha * alpha * sth2
	- 8. * a3 * g22 * g22 * alpha * alpha * sth2;

	F0 = -2. * beta * cth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * cth2
	+ a3 * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * sth2
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * sth2
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * cth2
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * sth2
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 2. * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * cth2
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * sth2
	+ a3 * intpow(g21, 3) * h21 * beta * beta * cth2
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * sth2);

	F1 = -2. * beta * cth * (-4. * a3 * g21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * beta * cth2 * alpha
	+ 2. * a3 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * h22 * g22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * g22 * g22 * beta * cth2 * alpha
	- 4. * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * g21 * g12 * h22 * alpha * sth2 * beta
	- 2. * a3 * g21 * g21 * g12 * h22 * beta * cth2 * alpha
	- 2. * a3 * g11 * h21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * g12 * g22 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * g12 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g11 * g21 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g11 * g21 * g22 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * g21 * g21 * beta * cth2 * alpha
	- 4. * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 2. * a3 * g11 * h21 * g22 * g22 * beta * cth2 * alpha
	+ 2. * a3 * g11 * h21 * alpha * sth2 * beta
	- 2. * a3 * g12 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * g21 * g21 * alpha * sth2 * beta
	+ 2. * a3 * g11 * h21 * g22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 2. * a3 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g22 * h22 * beta * cth2 * alpha)
	+ 4. * alpha * sth * (-a6 * a6 * h22 * g22
	- a6 * a6 * h21 * g21
	- a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * cth2
	+ a3 * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * sth2
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * sth2
	+ a3 * g21 * h21 * g22 * g22 * beta * beta * cth2
	- a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * sth2
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 2. * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * cth2
	+ a3 * g21 * g21 * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * sth2
	+ a3 * intpow(g21, 3) * h21 * beta * beta * cth2
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * sth2);

	F2 = -2. * beta * cth * (-(2 * a6 * a6 * h22 * g22)
	- (2 * a6 * a6 * h21 * g21)
	+ 2. * a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	- 2. * a3 * g11 * h21 * g12 * g22 * beta * beta * sth2
	+ 4. * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	- 2. * a3 * intpow(g22, 3) * h22 * beta * beta * cth2
	- 2. * a3 * g22 * h22 * beta * beta * cth2
	+ 2. * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 2. * a3 * g12 * h22 * beta * beta * sth2
	- 2. * a3 * g12 * g12 * h22 * g22 * beta * beta * sth2
	- 2. * a3 * g21 * h21 * (g22 * g22) * beta * beta * cth2
	+ 2. * a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g11 * h21 * beta * beta * cth * sth
	- 2. * a3 * g11 * h21 * beta * beta * sth2
	+ 2. * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g22 * h22 * beta * beta * cth * sth
	+ 4. * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * beta * beta * cth2
	- 2. * a3 * (g21 * g21) * g22 * h22 * beta * beta * cth2
	+ 2. * a3 * g21 * h21 * beta * beta * cth * sth
	+ 2. * a3 * g12 * h22 * beta * beta * cth * sth
	- 2. * a3 * g11 * g11 * h21 * g21 * beta * beta * sth2
	- 2. * a3 * intpow(g21, 3) * h21 * beta * beta * cth2
	- 2. * a3 * g11 * g21 * g12 * h22 * beta * beta * sth2
	+ 4. * a3 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a3 * (g21 * g21) * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * intpow(g21, 3) * h21 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 4. * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g12 * g12 * h22 * g22 * alpha * alpha * cth2
	+ 4. * a3 * (g21 * g21) * g12 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 4. * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * g21 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * g12 * g22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * (g22 * g22) * alpha * alpha * sth * cth
	+ 8. * a3 * g11 * h21 * (g21 * g21) * alpha * alpha * sth * cth
	+ 4. * a3 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * (g22 * g22) * alpha * alpha * sth2
	+ 4. * a3 * intpow(g22, 3) * h22 * alpha * alpha * sth2
	+ 4. * a3 * g11 * g11 * h21 * g21 * alpha * alpha * cth2
	+ 4. * a3 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 8. * a3 * g12 * h22 * (g22 * g22) * alpha * alpha * sth * cth)
	+ 4. * alpha * sth * (-4. * a3 * g21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * (g21 * g21) * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * beta * cth2 * alpha
	+ 2. * a3 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * h22 * (g22 * g22) * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * (g22 * g22) * beta * cth2 * alpha
	- 4. * a3 * g21 * h21 * (g22 * g22) * alpha * sth * beta * cth
	+ 2. * a3 * (g21 * g21) * g12 * h22 * alpha * sth2 * beta
	- 2. * a3 * (g21 * g21) * g12 * h22 * beta * cth2 * alpha
	- 2. * a3 * g11 * h21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * g12 * g22 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * g12 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g11 * g21 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g11 * g21 * g22 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * (g21 * g21) * beta * cth2 * alpha
	- 4. * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 2. * a3 * g11 * h21 * (g22 * g22) * beta * cth2 * alpha
	+ 2. * a3 * g11 * h21 * alpha * sth2 * beta
	- 2. * a3 * g12 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * (g21 * g21) * alpha * sth2 * beta
	+ 2. * a3 * g11 * h21 * (g22 * g22) * alpha * sth2 * beta
	+ 4. * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 2. * a3 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g22 * h22 * beta * cth2 * alpha)
	+ 2. * beta * cth * (-(a6 * a6 * h22 * g22)
	- (a6 * a6 * h21 * g21)
	- a3 * g11 * h21 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g11 * h21 * g12 * g22 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * (g21 * g21) * beta * beta * cth * sth
	+ a3 * intpow(g22, 3) * h22 * beta * beta * cth2
	+ a3 * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	+ a3 * g12 * h22 * beta * beta * sth2
	+ a3 * g12 * g12 * h22 * g22 * beta * beta * sth2
	+ a3 * g21 * h21 * (g22 * g22) * beta * beta * cth2
	- a3 * (g21 * g21) * g12 * h22 * beta * beta * cth * sth
	- a3 * g11 * h21 * beta * beta * cth * sth
	+ a3 * g11 * h21 * beta * beta * sth2
	- a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	- a3 * g22 * h22 * beta * beta * cth * sth
	- 2. * a3 * g12 * h22 * (g22 * g22) * beta * beta * cth * sth
	+ a3 * g21 * h21 * beta * beta * cth2
	+ a3 * (g21 * g21) * g22 * h22 * beta * beta * cth2
	- a3 * g21 * h21 * beta * beta * cth * sth
	- a3 * g12 * h22 * beta * beta * cth * sth
	+ a3 * g11 * g11 * h21 * g21 * beta * beta * sth2
	+ a3 * intpow(g21, 3) * h21 * beta * beta * cth2
	+ a3 * g11 * g21 * g12 * h22 * beta * beta * sth2);

	F3 = -2. * beta * cth * (4. * a3 * g21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	- 2. * a3 * g21 * h21 * alpha * sth2 * beta
	+ 2. * a3 * g21 * h21 * beta * cth2 * alpha
	- 2. * a3 * g12 * h22 * alpha * sth2 * beta
	- 4. * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g12 * h22 * g22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g12 * h22 * g22 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	- 2. * a3 * g21 * g21 * g12 * h22 * alpha * sth2 * beta
	+ 2. * a3 * g21 * g21 * g12 * h22 * beta * cth2 * alpha
	+ 2. * a3 * g11 * h21 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g22 * h22 * alpha * sth * beta * cth
	- 2. * a3 * g21 * h21 * g12 * g22 * alpha * sth2 * beta
	+ 2. * a3 * g21 * h21 * g12 * g22 * beta * cth2 * alpha
	- 4. * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	- 2. * a3 * g11 * g21 * g22 * h22 * alpha * sth2 * beta
	+ 2. * a3 * g11 * g21 * g22 * h22 * beta * cth2 * alpha
	- 4. * a3 * g12 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g21 * g21 * beta * cth2 * alpha
	+ 4. * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g11 * h21 * g22 * g22 * beta * cth2 * alpha
	- 2. * a3 * g11 * h21 * alpha * sth2 * beta
	+ 2. * a3 * g12 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * g21 * g21 * alpha * sth2 * beta
	- 2. * a3 * g11 * h21 * g22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	- 2. * a3 * g22 * h22 * alpha * sth2 * beta
	+ 2. * a3 * g22 * h22 * beta * cth2 * alpha)
	+ 4. * alpha * sth * (-2. * a6 * a6 * h22 * g22
	- 2. * a6 * a6 * h21 * g21
	+ 2. * a3 * g11 * h21 * g22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g11 * h21 * g12 * g22 * beta * beta * sth2
	+ 4. * a3 * g11 * h21 * g21 * g21 * beta * beta * cth * sth
	- 2. * a3 * intpow(g22, 3) * h22 * beta * beta * cth2
	- 2. * a3 * g22 * h22 * beta * beta * cth2
	+ 2. * a3 * g21 * h21 * g12 * g22 * beta * beta * cth * sth
	- 2. * a3 * g12 * h22 * beta * beta * sth2
	- 2. * a3 * g12 * g12 * h22 * g22 * beta * beta * sth2
	- 2. * a3 * g21 * h21 * g22 * g22 * beta * beta * cth2
	+ 2. * a3 * g21 * g21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g11 * h21 * beta * beta * cth * sth
	- 2. * a3 * g11 * h21 * beta * beta * sth2
	+ 2. * a3 * g11 * g21 * g22 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g22 * h22 * beta * beta * cth * sth
	+ 4. * a3 * g12 * h22 * g22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * beta * beta * cth2
	- 2. * a3 * g21 * g21 * g22 * h22 * beta * beta * cth2
	+ 2. * a3 * g21 * h21 * beta * beta * cth * sth
	+ 2. * a3 * g12 * h22 * beta * beta * cth * sth
	- 2. * a3 * g11 * g11 * h21 * g21 * beta * beta * sth2
	- 2. * a3 * intpow(g21, 3) * h21 * beta * beta * cth2
	- 2. * a3 * g11 * g21 * g12 * h22 * beta * beta * sth2
	+ 4. * a3 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a3 * g21 * g21 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * intpow(g21, 3) * h21 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * alpha * alpha * sth * cth
	+ 4. * a3 * g12 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g12 * g12 * h22 * g22 * alpha * alpha * cth2
	+ 4. * a3 * g21 * g21 * g12 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * alpha * alpha * sth * cth
	+ 4. * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * g21 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * g12 * g22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * h21 * g22 * g22 * alpha * alpha * sth * cth
	+ 8. * a3 * g11 * h21 * g21 * g21 * alpha * alpha * sth * cth
	+ 4. * a3 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * g22 * g22 * alpha * alpha * sth2
	+ 4. * a3 * intpow(g22, 3) * h22 * alpha * alpha * sth2
	+ 4. * a3 * g11 * g11 * h21 * g21 * alpha * alpha * cth2
	+ 4. * a3 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g11 * g21 * g22 * h22 * alpha * alpha * sth * cth
	+ 4. * a3 * g22 * h22 * alpha * alpha * sth * cth
	+ 8. * a3 * g12 * h22 * g22 * g22 * alpha * alpha * sth * cth)
	+ 2. * beta * cth * (-4. * a3 * g21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * intpow(g21, 3) * h21 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * beta * cth2 * alpha
	+ 2. * a3 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g12 * g12 * h22 * g22 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * h22 * g22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * g22 * g22 * beta * cth2 * alpha
	- 4. * a3 * g21 * h21 * g22 * g22 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * g21 * g12 * h22 * alpha * sth2 * beta
	- 2. * a3 * g21 * g21 * g12 * h22 * beta * cth2 * alpha
	- 2. * a3 * g11 * h21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g22 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g21 * h21 * g12 * g22 * alpha * sth2 * beta
	- 2. * a3 * g21 * h21 * g12 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g21 * g12 * h22 * alpha * sth * beta * cth
	+ 2. * a3 * g11 * g21 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g11 * g21 * g22 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g12 * g22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * g21 * g21 * beta * cth2 * alpha
	- 4. * a3 * intpow(g22, 3) * h22 * alpha * sth * beta * cth
	- 2. * a3 * g11 * h21 * g22 * g22 * beta * cth2 * alpha
	+ 2. * a3 * g11 * h21 * alpha * sth2 * beta
	- 2. * a3 * g12 * h22 * beta * cth2 * alpha
	+ 4. * a3 * g11 * h21 * g21 * g21 * alpha * sth2 * beta
	+ 2. * a3 * g11 * h21 * g22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * g11 * h21 * g21 * alpha * sth * beta * cth
	+ 2. * a3 * g22 * h22 * alpha * sth2 * beta
	- 2. * a3 * g22 * h22 * beta * cth2 * alpha);

	G0 = -beta * beta * cth2 * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 2. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	- 2. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	- 2. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G1 = -beta * beta * cth2 * (-4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta)
	+ 4. * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 2. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	- 2. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	- 2. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G2 = -beta * beta * cth2 * (-4. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * alpha * cth2
	- 2. * a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	+ 4. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a6 * a6 * h21 * h21
	+ 8. * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 2. * a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a6 * a6 * h22 * h22
	+ 8. * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 4. * a3 * g22 * g22 * h22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g21 * g21 * h21 * h21 * alpha * alpha * sth2
	- 4. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	- 2. * a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 8. * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 2. * a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	+ 4. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 8. * a3 * g21 * h21 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * alpha * cth2
	+ 8. * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	+ 4. * beta * cth * alpha * sth * (-4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta)
	- (-2. * beta * beta * cth2
	+ 4. * alpha * alpha * sth2) * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 2. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	- 2. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	- 2. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G3 = -beta * beta * cth2 * (4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	- 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta
	- 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha)
	+ 4. * beta * cth * alpha * sth * (-4. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * alpha * cth2
	- 2. * a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	+ 4. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a6 * a6 * h21 * h21
	+ 8. * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 2. * a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a6 * a6 * h22 * h22
	+ 8. * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 4. * a3 * g22 * g22 * h22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g21 * g21 * h21 * h21 * alpha * alpha * sth2
	- 4. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	- 2. * a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 8. * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 2. * a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	+ 4. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 8. * a3 * g21 * h21 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * alpha * cth2
	+ 8. * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- (-2. * beta * beta * cth2
	+ 4. * alpha * alpha * sth2) * (-4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta)
	- 4. * beta * cth * alpha * sth * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 2. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	- 2. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	- 2. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth);

	G4 = -2. * beta * beta * cth2 * (-a6 * a6 * h21 * h21
	- a6 * a6 * h22 * h22
	+ a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 2. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	- 2. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 2. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	+ a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	- 2. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth)
	+ 4. * beta * cth * alpha * sth * (4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	- 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	+ 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta
	- 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha)
	- (-2. * beta * beta * cth2
	+ 4. * alpha * alpha * sth2) * (-4. * a3 * g11 * h21 * g12 * h22 * beta * beta * sth2
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * alpha * cth2
	- 2. * a3 * g22 * g22 * h22 * h22 * beta * beta * cth2
	+ 4. * a3 * g12 * h22 * h22 * g22 * beta * beta * cth * sth
	- 2. * a6 * a6 * h21 * h21
	+ 8. * a3 * g21 * h21 * g12 * h22 * alpha * alpha * sth * cth
	- 2. * a3 * g11 * g11 * h21 * h21 * beta * beta * sth2
	- 2. * a6 * a6 * h22 * h22
	+ 8. * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sth * cth
	+ 4. * a3 * g11 * h21 * h21 * g21 * beta * beta * cth * sth
	+ 4. * a3 * g22 * g22 * h22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * alpha * cth2
	+ 4. * a3 * g21 * g21 * h21 * h21 * alpha * alpha * sth2
	- 4. * a3 * g21 * h21 * g22 * h22 * beta * beta * cth2
	- 2. * a3 * g12 * g12 * h22 * h22 * beta * beta * sth2
	+ 8. * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sth * cth
	- 2. * a3 * g21 * g21 * h21 * h21 * beta * beta * cth2
	+ 4. * a3 * g11 * h21 * g22 * h22 * beta * beta * cth * sth
	+ 8. * a3 * g21 * h21 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a3 * g21 * h21 * g12 * h22 * beta * beta * cth * sth
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * alpha * cth2
	+ 8. * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sth * cth)
	- 4. * beta * cth * alpha * sth * (-4. * a3 * g21 * h21 * g12 * h22 * beta * cth2 * alpha
	- 4. * a3 * g11 * h21 * h21 * g21 * beta * cth2 * alpha
	+ 4. * a3 * g11 * g11 * h21 * h21 * alpha * sth * beta * cth
	- 4. * a3 * g21 * g21 * h21 * h21 * alpha * sth * beta * cth
	+ 4. * a3 * g12 * g12 * h22 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g11 * h21 * g22 * h22 * alpha * sth2 * beta
	- 4. * a3 * g12 * h22 * h22 * g22 * beta * cth2 * alpha
	+ 4. * a3 * g12 * h22 * h22 * g22 * alpha * sth2 * beta
	- 4. * a3 * g11 * h21 * g22 * h22 * beta * cth2 * alpha
	- 8. * a3 * g21 * h21 * g22 * h22 * alpha * sth * beta * cth
	- 4. * a3 * g22 * g22 * h22 * h22 * alpha * sth * beta * cth
	+ 8. * a3 * g11 * h21 * g12 * h22 * alpha * sth * beta * cth
	+ 4. * a3 * g21 * h21 * g12 * h22 * alpha * sth2 * beta
	+ 4. * a3 * g11 * h21 * h21 * g21 * alpha * sth2 * beta);

	H0 = -beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	- 2. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * sth2);

	H1 = -beta * beta * sth * cth * (4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta
	- 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta)
	- (-2. * beta * sth2 * alpha
	+ 2. * alpha * cth2 * beta) * (a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	- 2. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * sth2);

	H2 = -beta * beta * sth * cth * (-2. * a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- 2. * a4 * a6 * h11 * h21
	- 2. * a4 * a6 * h12 * h22
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 2. * a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g12 * g12 * h12 * h22 * beta * beta * sth2
	+ 4. * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g21 * g21 * h11 * h21 * alpha * alpha * sth2
	+ 8. * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 8. * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g22 * g22 * h12 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * g11 * h11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-2. * beta * sth2 * alpha
	+ 2. * alpha * cth2 * beta) * (4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta
	- 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta)
	- (-2. * beta * beta * sth * cth
	- 4. * alpha * alpha * cth * sth) * (a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	- 2. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * sth2);

	H3 = -beta * beta * sth * cth * (4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	+ 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta)
	- (-2. * beta * sth2 * alpha
	+ 2. * alpha * cth2 * beta) * (-2. * a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- 2. * a4 * a6 * h11 * h21
	- 2. * a4 * a6 * h12 * h22
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 2. * a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g12 * g12 * h12 * h22 * beta * beta * sth2
	+ 4. * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g21 * g21 * h11 * h21 * alpha * alpha * sth2
	+ 8. * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 8. * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g22 * g22 * h12 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * g11 * h11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-2. * beta * beta * sth * cth
	- 4. * alpha * alpha * cth * sth) * (4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta
	- 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta)
	- (-2. * alpha * cth2 * beta
	+ 2. * beta * sth2 * alpha) * (a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	- 2. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * sth2);

	H4 = -2. * beta * beta * sth * cth * (a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- a4 * a6 * h11 * h21
	- a4 * a6 * h12 * h22
	- a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	+ a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	- 2. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	- a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	+ a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	- a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * h22 * beta * beta * sth2)
	- (-2. * beta * sth2 * alpha
	+ 2. * alpha * cth2 * beta) * (4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	- 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	+ 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta)
	- (-2. * beta * beta * sth * cth
	- 4. * alpha * alpha * cth * sth) * (-2. * a5 * g11 * g11 * h11 * h21 * beta * beta * sth2
	- 2. * a4 * a6 * h11 * h21
	- 2. * a4 * a6 * h12 * h22
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * beta * sth * cth
	- 2. * a5 * g22 * h12 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g22 * g22 * h12 * h22 * beta * beta * cth2
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g12 * h22 * beta * beta * sth2
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g21 * g21 * h11 * h21 * beta * beta * cth2
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * g21 * h11 * g22 * h22 * beta * beta * cth2
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * beta * sth * cth
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g12 * g12 * h12 * h22 * beta * beta * sth2
	+ 4. * a5 * g11 * h11 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g21 * g21 * h11 * h21 * alpha * alpha * sth2
	+ 8. * a5 * g11 * h11 * g21 * h21 * alpha * alpha * cth * sth
	+ 8. * a5 * g12 * h12 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g22 * g22 * h12 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h12 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * g11 * h11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h12 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g21 * h11 * g12 * h22 * alpha * alpha * cth * sth)
	- (-2. * alpha * cth2 * beta
	+ 2. * beta * sth2 * alpha) * (4. * a5 * g11 * g11 * h11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g12 * h12 * g11 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * h12 * g21 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * h11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * g11 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * h11 * g22 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * h21 * alpha * cth2 * beta
	- 4. * a5 * g22 * g22 * h12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g22 * h22 * alpha * cth2 * beta);

	J0 = -beta * cth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * sth2
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * cth2
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * sth2
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12_3 * h22 * beta * beta * sth2
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * cth2
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * cth2
	+ a5 * g11_3 * h21 * beta * beta * sth2
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * cth2
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * cth2
	+ a5 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J1 = -beta * cth * (4. * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * g11 * g22 * h22 * alpha * cth2 * beta
	- 2. * a5 * g21 * g11 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * g11 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g11 * g11 * g21 * h21 * alpha * cth2 * beta
	+ 4. * a5 * g11 * g11 * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12_3 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g11_3 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g12 * g12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * g12 * g11 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * g12 * g11 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h21 * alpha * cth2 * beta
	- 2. * a5 * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g12 * g12 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * g12 * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g12 * g12 * g22 * h22 * beta * sth2 * alpha
	+ 2. * a5 * g11 * h21 * beta * sth2 * alpha
	+ 2. * a5 * g11 * g11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 2. * alpha * sth * (-a4 * a6 * g12 * h22
	- a4 * a6 * g11 * h21
	+ a5 * g12 * h22 * beta * beta * sth2
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * cth2
	+ a5 * g12 * g12 * g11 * h21 * beta * beta * sth2
	- a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12_3 * h22 * beta * beta * sth2
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * cth2
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * cth2
	+ a5 * g11_3 * h21 * beta * beta * sth2
	+ a5 * g11 * g11 * g12 * h22 * beta * beta * sth2
	- a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * cth2
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * cth2
	+ a5 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth);

	J2 = -beta * cth * (-(2 * a4 * a6 * g12 * h22)
	- (2 * a4 * a6 * g11 * h21)
	- 2. * a5 * g12 * h22 * beta * beta * sth2
	- 2. * a5 * g22 * g12 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * (g12 * g12) * g11 * h21 * beta * beta * sth2
	+ 2. * a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g11 * h21 * beta * beta * sth2
	+ 4. * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g12 * h22 * beta * beta * sth * cth
	+ 2. * a5 * g11 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 2. * a5 * g12_3 * h22 * beta * beta * sth2
	- 2. * a5 * g22 * g22 * g12 * h22 * beta * beta * cth2
	+ 2. * a5 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g11_3 * h21 * beta * beta * sth2
	- 2. * a5 * (g11 * g11) * g12 * h22 * beta * beta * sth2
	+ 2. * a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g21 * g21 * g11 * h21 * beta * beta * cth2
	- 2. * a5 * g21 * g11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g22 * h22 * beta * beta * cth2
	+ 4. * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth
	+ 4. * a5 * (g11 * g11) * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * (g11 * g11) * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 8. * a5 * (g11 * g11) * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12_3 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g22 * g22 * g12 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11_3 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * g12 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g11 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * (g12 * g12) * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * (g12 * g12) * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g22 * h22 * alpha * alpha * sth2
	+ 8. * a5 * (g12 * g12) * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g21 * g11 * h21 * alpha * alpha * sth2)
	+ 2. * alpha * sth * (4. * a5 * (g11 * g11) * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * (g11 * g11) * g22 * h22 * alpha * cth2 * beta
	- 2. * a5 * g21 * g11 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * g11 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * (g11 * g11) * g21 * h21 * alpha * cth2 * beta
	+ 4. * a5 * (g11 * g11) * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12_3 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g11_3 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * (g12 * g12) * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * g12 * g11 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * g12 * g11 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h21 * alpha * cth2 * beta
	- 2. * a5 * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * (g12 * g12) * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * (g12 * g12) * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * (g12 * g12) * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * (g12 * g12) * g22 * h22 * beta * sth2 * alpha
	+ 2. * a5 * g11 * h21 * beta * sth2 * alpha
	+ 2. * a5 * (g11 * g11) * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ beta * cth * (-(a4 * a6 * g12 * h22)
	- (a4 * a6 * g11 * h21)
	+ a5 * g12 * h22 * beta * beta * sth2
	+ a5 * g22 * g12 * g21 * h21 * beta * beta * cth2
	+ a5 * (g12 * g12) * g11 * h21 * beta * beta * sth2
	- a5 * (g12 * g12) * g21 * h21 * beta * beta * sth * cth
	- a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	- a5 * g22 * h22 * beta * beta * sth * cth
	+ a5 * g11 * h21 * beta * beta * sth2
	- 2. * a5 * (g11 * g11) * g21 * h21 * beta * beta * sth * cth
	- a5 * g12 * h22 * beta * beta * sth * cth
	- a5 * g11 * h21 * beta * beta * sth * cth
	- a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	+ a5 * g12_3 * h22 * beta * beta * sth2
	+ a5 * g22 * g22 * g12 * h22 * beta * beta * cth2
	- a5 * g21 * h21 * beta * beta * sth * cth
	+ a5 * g21 * h21 * beta * beta * cth2
	+ a5 * g11_3 * h21 * beta * beta * sth2
	+ a5 * (g11 * g11) * g12 * h22 * beta * beta * sth2
	- a5 * (g11 * g11) * g22 * h22 * beta * beta * sth * cth
	+ a5 * g21 * g21 * g11 * h21 * beta * beta * cth2
	+ a5 * g21 * g11 * g22 * h22 * beta * beta * cth2
	+ a5 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * (g12 * g12) * g22 * h22 * beta * beta * sth * cth);

	J3 = -beta * cth * (-4. * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * g11 * g22 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * g11 * g12 * h22 * alpha * cth2 * beta
	- 2. * a5 * g21 * g11 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g11 * g11 * g21 * h21 * alpha * cth2 * beta
	- 4. * a5 * g11 * g11 * g21 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	- 4. * a5 * g12_3 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g12 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g21 * h21 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g11_3 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * h22 * alpha * cth2 * beta
	- 2. * a5 * g12 * g12 * g21 * h21 * beta * sth2 * alpha
	+ 2. * a5 * g22 * g12 * g11 * h21 * alpha * cth2 * beta
	- 2. * a5 * g22 * g12 * g11 * h21 * beta * sth2 * alpha
	- 4. * a5 * g11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * g12 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g22 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h21 * alpha * cth2 * beta
	- 2. * a5 * g21 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g22 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * g22 * h22 * alpha * cth2 * beta
	- 4. * a5 * g12 * g12 * g22 * h22 * beta * sth2 * alpha
	- 2. * a5 * g11 * h21 * beta * sth2 * alpha
	- 2. * a5 * g11 * g11 * g22 * h22 * beta * sth2 * alpha
	+ 4. * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth)
	+ 2. * alpha * sth * (-2. * a4 * a6 * g12 * h22
	- 2. * a4 * a6 * g11 * h21
	- 2. * a5 * g12 * h22 * beta * beta * sth2
	- 2. * a5 * g22 * g12 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g12 * g12 * g11 * h21 * beta * beta * sth2
	+ 2. * a5 * g12 * g12 * g21 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g22 * g12 * g11 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g11 * h21 * beta * beta * sth2
	+ 4. * a5 * g11 * g11 * g21 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g12 * h22 * beta * beta * sth * cth
	+ 2. * a5 * g11 * h21 * beta * beta * sth * cth
	+ 2. * a5 * g21 * g11 * g12 * h22 * beta * beta * sth * cth
	- 2. * a5 * g12_3 * h22 * beta * beta * sth2
	- 2. * a5 * g22 * g22 * g12 * h22 * beta * beta * cth2
	+ 2. * a5 * g21 * h21 * beta * beta * sth * cth
	- 2. * a5 * g21 * h21 * beta * beta * cth2
	- 2. * a5 * g11_3 * h21 * beta * beta * sth2
	- 2. * a5 * g11 * g11 * g12 * h22 * beta * beta * sth2
	+ 2. * a5 * g11 * g11 * g22 * h22 * beta * beta * sth * cth
	- 2. * a5 * g21 * g21 * g11 * h21 * beta * beta * cth2
	- 2. * a5 * g21 * g11 * g22 * h22 * beta * beta * cth2
	- 2. * a5 * g22 * h22 * beta * beta * cth2
	+ 4. * a5 * g12 * g12 * g22 * h22 * beta * beta * sth * cth
	+ 4. * a5 * g11 * g11 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * g11 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g11 * g12 * h22 * alpha * alpha * cth * sth
	+ 8. * a5 * g11 * g11 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12_3 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g22 * g22 * g12 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h22 * alpha * alpha * cth2
	+ 4. * a5 * g11_3 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * g12 * g21 * h21 * alpha * alpha * sth2
	+ 4. * a5 * g22 * g12 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g11 * g22 * h22 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * g12 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h21 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * g12 * g11 * h21 * alpha * alpha * cth2
	+ 4. * a5 * g22 * h22 * alpha * alpha * sth2
	+ 8. * a5 * g12 * g12 * g22 * h22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * g21 * g11 * h21 * alpha * alpha * sth2)
	+ beta * cth * (4. * a5 * g11 * g11 * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g11 * g11 * g22 * h22 * alpha * cth2 * beta
	- 2. * a5 * g21 * g11 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * g11 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g11 * g11 * g21 * h21 * alpha * cth2 * beta
	+ 4. * a5 * g11 * g11 * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g22 * g12 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12_3 * h22 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * h21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h22 * alpha * cth * beta * sth
	+ 4. * a5 * g11_3 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g12 * h22 * alpha * cth2 * beta
	+ 2. * a5 * g12 * g12 * g21 * h21 * beta * sth2 * alpha
	- 2. * a5 * g22 * g12 * g11 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * g12 * g11 * h21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h21 * alpha * cth2 * beta
	- 2. * a5 * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g12 * g12 * g11 * h21 * alpha * cth * beta * sth
	- 2. * a5 * g12 * g12 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g22 * g22 * g12 * h22 * alpha * cth * beta * sth
	- 2. * a5 * g21 * h21 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h21 * beta * sth2 * alpha
	- 4. * a5 * g21 * g11 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g22 * h22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * g12 * g22 * h22 * alpha * cth2 * beta
	+ 4. * a5 * g12 * g12 * g22 * h22 * beta * sth2 * alpha
	+ 2. * a5 * g11 * h21 * beta * sth2 * alpha
	+ 2. * a5 * g11 * g11 * g22 * h22 * beta * sth2 * alpha
	- 4. * a5 * g21 * g21 * g11 * h21 * alpha * cth * beta * sth);

	K0 = -beta * sth * (a5 * intpow(g21, 3) * h11 * beta * beta * cth2
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * sth2
	+ a5 * g21 * h11 * beta * beta * cth2
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * sth2
	+ a5 * g22 * h12 * beta * beta * cth2
	+ a5 * intpow(g22, 3) * h12 * beta * beta * cth2
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * sth2
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * cth2
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * cth2
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * sth2
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * sth2
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K1 = -beta * sth * (2. * a5 * g22 * h12 * g11 * g21 * beta * sth2 * alpha
	+ 2. * a5 * g12 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * g21 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * g22 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * beta * sth2 * alpha
	- 4. * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * g21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g21 * g21 * beta * sth2 * alpha
	- 4. * a5 * g11 * h11 * g21 * g21 * alpha * cth2 * beta
	+ 2. * a5 * g11 * h11 * g22 * g22 * beta * sth2 * alpha
	- 2. * a5 * g11 * h11 * g22 * g22 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * g22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * alpha * cth2 * beta
	- 2. * a5 * g12 * h12 * g21 * g21 * alpha * cth2 * beta
	- 4. * a5 * g12 * h12 * g22 * g22 * alpha * cth2 * beta
	- 4. * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * beta * sth2 * alpha
	- 2. * a5 * g21 * h11 * alpha * cth2 * beta
	+ 4. * a5 * g12 * h12 * g22 * g22 * beta * sth2 * alpha)
	- 2. * alpha * cth * (a5 * intpow(g21, 3) * h11 * beta * beta * cth2
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * sth2
	+ a5 * g21 * h11 * beta * beta * cth2
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * sth2
	+ a5 * g22 * h12 * beta * beta * cth2
	+ a5 * intpow(g22, 3) * h12 * beta * beta * cth2
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * sth2
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * cth2
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * cth2
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * sth2
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * sth2
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K2 = -beta * sth * (-2. * a5 * intpow(g21, 3) * h11 * beta * beta * cth2
	+ 2. * a5 * g12 * h12 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * beta * beta * sth2
	- 2. * a5 * g21 * h11 * beta * beta * cth2
	+ 2. * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 2. * a5 * g11 * h11 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * beta * beta * sth2
	- 2. * a5 * g22 * h12 * beta * beta * cth2
	- 2. * a5 * intpow(g22, 3) * h12 * beta * beta * cth2
	- 2. * a5 * g11 * g11 * h11 * g21 * beta * beta * sth2
	+ 2. * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 2. * a5 * g22 * h12 * g21 * g21 * beta * beta * cth2
	- 2. * a5 * g21 * h11 * g22 * g22 * beta * beta * cth2
	+ 2. * a5 * g21 * h11 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g11 * g21 * beta * beta * sth2
	+ 2. * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 4. * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g12 * g22 * beta * beta * sth2
	+ 2. * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 2. * a5 * g12 * g12 * h12 * g22 * beta * beta * sth2
	+ 2. * a5 * g22 * h12 * beta * beta * sth * cth
	+ 4. * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 2. * a4 * a6 * h11 * g21
	- 2. * a4 * a6 * h12 * g22
	+ 4. * a5 * intpow(g22, 3) * h12 * alpha * alpha * sth2
	+ 4. * a5 * intpow(g21, 3) * h11 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * g12 * g22 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * alpha * alpha * sth2
	+ 4. * a5 * g12 * g12 * h12 * g22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * h12 * g11 * g21 * alpha * alpha * cth2
	+ 8. * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * alpha * alpha * cth2
	+ 4. * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * g21 * g21 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * g11 * h11 * g21 * alpha * alpha * cth2
	+ 8. * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g22 * g22 * alpha * alpha * sth2)
	- 2. * alpha * cth * (2. * a5 * g22 * h12 * g11 * g21 * beta * sth2 * alpha
	+ 2. * a5 * g12 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * g21 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * g22 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * beta * sth2 * alpha
	- 4. * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * g21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g21 * g21 * beta * sth2 * alpha
	- 4. * a5 * g11 * h11 * g21 * g21 * alpha * cth2 * beta
	+ 2. * a5 * g11 * h11 * g22 * g22 * beta * sth2 * alpha
	- 2. * a5 * g11 * h11 * g22 * g22 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * g22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * alpha * cth2 * beta
	- 2. * a5 * g12 * h12 * g21 * g21 * alpha * cth2 * beta
	- 4. * a5 * g12 * h12 * g22 * g22 * alpha * cth2 * beta
	- 4. * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * beta * sth2 * alpha
	- 2. * a5 * g21 * h11 * alpha * cth2 * beta
	+ 4. * a5 * g12 * h12 * g22 * g22 * beta * sth2 * alpha)
	+ beta * sth * (a5 * intpow(g21, 3) * h11 * beta * beta * cth2
	- a5 * g12 * h12 * beta * beta * sth * cth
	+ a5 * g11 * h11 * beta * beta * sth2
	+ a5 * g21 * h11 * beta * beta * cth2
	- a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	- a5 * g11 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * beta * beta * sth2
	+ a5 * g22 * h12 * beta * beta * cth2
	+ a5 * intpow(g22, 3) * h12 * beta * beta * cth2
	+ a5 * g11 * g11 * h11 * g21 * beta * beta * sth2
	- a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	+ a5 * g22 * h12 * g21 * g21 * beta * beta * cth2
	+ a5 * g21 * h11 * g22 * g22 * beta * beta * cth2
	- a5 * g21 * h11 * beta * beta * sth * cth
	+ a5 * g12 * h12 * g11 * g21 * beta * beta * sth2
	- a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * h11 * g12 * g22 * beta * beta * sth2
	- a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * h12 * g22 * beta * beta * sth2
	- a5 * g22 * h12 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- a4 * a6 * h11 * g21
	- a4 * a6 * h12 * g22);

	K3 = -beta * sth * (-2. * a5 * g22 * h12 * g11 * g21 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * beta * sth2 * alpha
	+ 2. * a5 * g22 * h12 * g11 * g21 * alpha * cth2 * beta
	+ 4. * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	- 4. * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	+ 4. * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	- 4. * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	- 4. * a5 * g12 * h12 * alpha * cth * beta * sth
	+ 4. * a5 * g21 * h11 * alpha * cth * beta * sth
	- 2. * a5 * g21 * h11 * g12 * g22 * beta * sth2 * alpha
	+ 2. * a5 * g12 * h12 * alpha * cth2 * beta
	+ 4. * a5 * g22 * h12 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * beta * sth2 * alpha
	+ 4. * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * alpha * cth2 * beta
	- 2. * a5 * g12 * h12 * g21 * g21 * beta * sth2 * alpha
	- 4. * a5 * g11 * h11 * alpha * cth * beta * sth
	- 4. * a5 * g11 * h11 * g21 * g21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h11 * g21 * g21 * alpha * cth2 * beta
	- 2. * a5 * g11 * h11 * g22 * g22 * beta * sth2 * alpha
	+ 2. * a5 * g11 * h11 * g22 * g22 * alpha * cth2 * beta
	+ 2. * a5 * g21 * h11 * g12 * g22 * alpha * cth2 * beta
	- 4. * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	- 2. * a5 * g22 * h12 * beta * sth2 * alpha
	+ 2. * a5 * g22 * h12 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * g21 * alpha * cth2 * beta
	+ 4. * a5 * g12 * h12 * g22 * g22 * alpha * cth2 * beta
	+ 4. * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	- 2. * a5 * g21 * h11 * beta * sth2 * alpha
	+ 2. * a5 * g21 * h11 * alpha * cth2 * beta
	- 4. * a5 * g12 * h12 * g22 * g22 * beta * sth2 * alpha)
	- 2. * alpha * cth * (-2. * a5 * intpow(g21, 3) * h11 * beta * beta * cth2
	+ 2. * a5 * g12 * h12 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * beta * beta * sth2
	- 2. * a5 * g21 * h11 * beta * beta * cth2
	+ 2. * a5 * g22 * h12 * g11 * g21 * beta * beta * sth * cth
	+ 2. * a5 * g11 * h11 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * beta * beta * sth2
	- 2. * a5 * g22 * h12 * beta * beta * cth2
	- 2. * a5 * intpow(g22, 3) * h12 * beta * beta * cth2
	- 2. * a5 * g11 * g11 * h11 * g21 * beta * beta * sth2
	+ 2. * a5 * g21 * h11 * g12 * g22 * beta * beta * sth * cth
	- 2. * a5 * g22 * h12 * g21 * g21 * beta * beta * cth2
	- 2. * a5 * g21 * h11 * g22 * g22 * beta * beta * cth2
	+ 2. * a5 * g21 * h11 * beta * beta * sth * cth
	- 2. * a5 * g12 * h12 * g11 * g21 * beta * beta * sth2
	+ 2. * a5 * g12 * h12 * g21 * g21 * beta * beta * sth * cth
	+ 4. * a5 * g11 * h11 * g21 * g21 * beta * beta * sth * cth
	- 2. * a5 * g11 * h11 * g12 * g22 * beta * beta * sth2
	+ 2. * a5 * g11 * h11 * g22 * g22 * beta * beta * sth * cth
	- 2. * a5 * g12 * g12 * h12 * g22 * beta * beta * sth2
	+ 2. * a5 * g22 * h12 * beta * beta * sth * cth
	+ 4. * a5 * g12 * h12 * g22 * g22 * beta * beta * sth * cth
	- 2. * a4 * a6 * h11 * g21
	- 2. * a4 * a6 * h12 * g22
	+ 4. * a5 * intpow(g22, 3) * h12 * alpha * alpha * sth2
	+ 4. * a5 * intpow(g21, 3) * h11 * alpha * alpha * sth2
	+ 4. * a5 * g22 * h12 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * alpha * alpha * cth2
	+ 4. * a5 * g22 * h12 * g11 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * g12 * g22 * alpha * alpha * cth2
	+ 4. * a5 * g12 * h12 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * alpha * alpha * sth2
	+ 4. * a5 * g12 * g12 * h12 * g22 * alpha * alpha * cth2
	+ 4. * a5 * g11 * h11 * alpha * alpha * cth * sth
	+ 4. * a5 * g12 * h12 * g11 * g21 * alpha * alpha * cth2
	+ 8. * a5 * g11 * h11 * g21 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * h11 * alpha * alpha * cth2
	+ 4. * a5 * g11 * h11 * g22 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g12 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * alpha * alpha * cth * sth
	+ 4. * a5 * g22 * h12 * g21 * g21 * alpha * alpha * sth2
	+ 4. * a5 * g12 * h12 * g21 * g21 * alpha * alpha * cth * sth
	+ 4. * a5 * g11 * g11 * h11 * g21 * alpha * alpha * cth2
	+ 8. * a5 * g12 * h12 * g22 * g22 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * alpha * alpha * cth * sth
	+ 4. * a5 * g21 * h11 * g22 * g22 * alpha * alpha * sth2)
	+ beta * sth * (2. * a5 * g22 * h12 * g11 * g21 * beta * sth2 * alpha
	+ 2. * a5 * g12 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * g11 * g21 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * g21 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * g12 * h12 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * g22 * g22 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * g11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * g11 * h11 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g12 * h12 * alpha * cth * beta * sth
	- 4. * a5 * g21 * h11 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * g12 * g22 * beta * sth2 * alpha
	- 2. * a5 * g12 * h12 * alpha * cth2 * beta
	- 4. * a5 * g22 * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * h11 * beta * sth2 * alpha
	- 4. * a5 * intpow(g21, 3) * h11 * alpha * cth * beta * sth
	- 2. * a5 * g11 * h11 * alpha * cth2 * beta
	+ 2. * a5 * g12 * h12 * g21 * g21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * h11 * alpha * cth * beta * sth
	+ 4. * a5 * g11 * h11 * g21 * g21 * beta * sth2 * alpha
	- 4. * a5 * g11 * h11 * g21 * g21 * alpha * cth2 * beta
	+ 2. * a5 * g11 * h11 * g22 * g22 * beta * sth2 * alpha
	- 2. * a5 * g11 * h11 * g22 * g22 * alpha * cth2 * beta
	- 2. * a5 * g21 * h11 * g12 * g22 * alpha * cth2 * beta
	+ 4. * a5 * g11 * h11 * g12 * g22 * alpha * cth * beta * sth
	+ 2. * a5 * g22 * h12 * beta * sth2 * alpha
	- 2. * a5 * g22 * h12 * alpha * cth2 * beta
	- 2. * a5 * g12 * h12 * g21 * g21 * alpha * cth2 * beta
	- 4. * a5 * g12 * h12 * g22 * g22 * alpha * cth2 * beta
	- 4. * a5 * intpow(g22, 3) * h12 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * h11 * beta * sth2 * alpha
	- 2. * a5 * g21 * h11 * alpha * cth2 * beta
	+ 4. * a5 * g12 * h12 * g22 * g22 * beta * sth2 * alpha);

L0 = a4 * a6 * g11 * g21
	+ a4 * a6 * g12 * g22
	+ a4 * a6
	- a5 * g22 * g22 * beta * beta * cth2
	- a5 * g21 * g21 * beta * beta * cth2
	+ a5 * g11 * g11 * g22 * g22 * beta * beta * sth * cth
	- a5 * beta * beta * sth2
	- a5 * beta * beta * cth2
	- a5 * g12_3 * g22 * beta * beta * sth2
	- a5 * g12 * g12 * g11 * g21 * beta * beta * sth2
	- a5 * intpow(g21, 3) * g11 * beta * beta * cth2
	- a5 * g22 * g12 * beta * beta * cth2
	- a5 * g21 * g11 * beta * beta * sth2
	- a5 * intpow(g22, 3) * g12 * beta * beta * cth2
	+ 2. * a5 * g21 * g11 * beta * beta * sth * cth
	- a5 * g22 * g12 * beta * beta * sth2
	+ 2. * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * beta * beta * cth2
	+ a5 * g12 * g12 * beta * beta * sth * cth
	- a5 * g12 * g12 * beta * beta * sth2
	+ 2. * a5 * g11 * g11 * g21 * g21 * beta * beta * sth * cth
	+ 2. * a5 * g22 * g12 * beta * beta * sth * cth
	+ 2. * a5 * beta * beta * sth * cth
	- a5 * g22 * g12 * g21 * g21 * beta * beta * cth2
	+ 2. * a5 * g12 * g12 * g22 * g22 * beta * beta * sth * cth
	+ a5 * g12 * g12 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g11 * g11 * beta * beta * sth * cth
	- a5 * g11 * g11 * g12 * g22 * beta * beta * sth2
	- a5 * g11_3 * g21 * beta * beta * sth2
	+ a5 * g21 * g21 * beta * beta * sth * cth
	+ a5 * g22 * g22 * beta * beta * sth * cth
	- a5 * g21 * g11 * g22 * g22 * beta * beta * cth2
	- a5 * g11 * g11 * beta * beta * sth2;

L1 = 4. * a5 * intpow(g22, 3) * g12 * alpha * cth * beta * sth
	+ 4. * a5 * intpow(g21, 3) * g11 * alpha * cth * beta * sth
	+ 4. * a5 * g21 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g21 * g11 * g12 * g22 * alpha * cth2 * beta
	- 4. * a5 * g21 * g11 * g12 * g22 * beta * sth2 * alpha
	+ 4. * a5 * g21 * g11 * alpha * cth2 * beta
	- 4. * a5 * g21 * g11 * beta * sth2 * alpha
	- 4. * a5 * g12 * g12 * g11 * g21 * alpha * cth * beta * sth
	- 4. * a5 * g12_3 * g22 * alpha * cth * beta * sth
	+ 2. * a5 * g11 * g11 * alpha * cth2 * beta
	- 2. * a5 * g11 * g11 * beta * sth2 * alpha
	+ 4. * a5 * g12 * g12 * g22 * g22 * alpha * cth2 * beta
	+ 4. * a5 * alpha * cth2 * beta
	- 4. * a5 * beta * sth2 * alpha
	- 4. * a5 * g12 * g12 * alpha * cth * beta * sth
	+ 4. * a5 * g22 * g12 * g21 * g21 * alpha * cth * beta * sth
	+ 4. * a5 * g22 * g12 * alpha * cth2 * beta
	- 4. * a5 * g22 * g12 * beta * sth2 * alpha
	- 2. * a5 * g21 * g21 * beta * sth2 * alpha
	- 4. * a5 * g11 * g11 * g21 * g21 * beta * sth2 * alpha
	+ 4. * a5 * g11 * g11 * g21 * g21 * alpha * cth2 * beta
	+ 4. * a5 * g22 * g22 * alpha * cth * beta * sth
	+ 2. * a5 * g12 * g12 * alpha * cth2 * beta
	- 2. * a5 * g12 * g12 * beta * sth2 * alpha
	+ 2. * a5 * g11 * g11 * g22 * g22 * alpha * cth2 * beta
	- 2. * a5 * g11 * g11 * g22 * g22 * beta * sth2 * alpha
	+ 4. * a5 * g21 * g11 * g22 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g11_3 * g21 * alpha * cth * beta * sth
	+ 2. * a5 * g21 * g21 * alpha * cth2 * beta
	+ 2. * a5 * g22 * g22 * alpha * cth2 * beta
	- 2. * a5 * g22 * g22 * beta * sth2 * alpha
	- 4. * a5 * g11 * g11 * g12 * g22 * alpha * cth * beta * sth
	- 4. * a5 * g11 * g11 * alpha * cth * beta * sth
	- 4. * a5 * g12 * g12 * g22 * g22 * beta * sth2 * alpha
	+ 2. * a5 * g12 * g12 * g21 * g21 * alpha * cth2 * beta
	- 2. * a5 * g12 * g12 * g21 * g21 * beta * sth2 * alpha;

L2 = (2 * a4 * a6 * g11 * g21)
	+ (2 * a4 * a6 * g12 * g22)
	+ (2 * a4 * a6)
	+ 2. * a5 * (g22 * g22) * beta * beta * cth2
	+ 2. * a5 * (g21 * g21) * beta * beta * cth2
	- 2. * a5 * (g11 * g11) * (g22 * g22) * beta * beta * sth * cth
	+ 2. * a5 * beta * beta * sth2
	+ 2. * a5 * beta * beta * cth2
	+ 2. * a5 * g12_3 * g22 * beta * beta * sth2
	+ 2. * a5 * (g12 * g12) * g11 * g21 * beta * beta * sth2
	+ 2. * a5 * intpow(g21, 3) * g11 * beta * beta * cth2
	+ 2. * a5 * g22 * g12 * beta * beta * cth2
	+ 2. * a5 * g21 * g11 * beta * beta * sth2
	+ 2. * a5 * intpow(g22, 3) * g12 * beta * beta * cth2
	- 4. * a5 * g21 * g11 * beta * beta * sth * cth
	+ 2. * a5 * g22 * g12 * beta * beta * sth2
	- 4. * a5 * g21 * g11 * g12 * g22 * beta * beta * sth * cth
	+ 2. * a5 * g21 * g11 * beta * beta * cth2
	- 2. * a5 * (g12 * g12) * beta * beta * sth * cth
	+ 2. * a5 * (g12 * g12) * beta * beta * sth2
	- 4. * a5 * (g11 * g11) * (g21 * g21) * beta * beta * sth * cth
	- 4. * a5 * g22 * g12 * beta * beta * sth * cth
	- 4. * a5 * beta * beta * sth * cth
	+ 2. * a5 * g22 * g12 * (g21 * g21) * beta * beta * cth2
	- 4. * a5 * (g12 * g12) * (g22 * g22) * beta * beta * sth * cth
	- 2. * a5 * (g12 * g12) * (g21 * g21) * beta * beta * sth * cth
	- 2. * a5 * (g11 * g11) * beta * beta * sth * cth
	+ 2. * a5 * (g11 * g11) * g12 * g22 * beta * beta * sth2
	+ 2. * a5 * g11_3 * g21 * beta * beta * sth2
	- 2. * a5 * (g21 * g21) * beta * beta * sth * cth
	- 2. * a5 * (g22 * g22) * beta * beta * sth * cth
	+ 2. * a5 * g21 * g11 * (g22 * g22) * beta * beta * cth2
	+ 2. * a5 * (g11 * g11) * beta * beta * sth2
	- 4. * a5 * alpha * alpha * cth2
	- 4. * a5 * alpha * alpha * sth2
	- 4. * a5 * intpow(g22, 3) * g12 * alpha * alpha * sth2
	- 4. * a5 * (g21 * g21) * alpha * alpha * sth2
	- 8. * a5 * (g11 * g11) * (g21 * g21) * alpha * alpha * cth * sth
	- 4. * a5 * g22 * g12 * alpha * alpha * cth2
	- 4. * a5 * (g11 * g11) * alpha * alpha * cth2
	- 4. * a5 * (g22 * g22) * alpha * alpha * sth2
	- 8. * a5 * alpha * alpha * cth * sth
	- 4. * a5 * intpow(g21, 3) * g11 * alpha * alpha * sth2
	- 4. * a5 * g22 * g12 * alpha * alpha * sth2
	- 4. * a5 * (g12 * g12) * g11 * g21 * alpha * alpha * cth2
	- 4. * a5 * g21 * g11 * alpha * alpha * cth2
	- 4. * a5 * g12_3 * g22 * alpha * alpha * cth2
	- 4. * a5 * (g11 * g11) * alpha * alpha * cth * sth
	- 8. * a5 * (g12 * g12) * (g22 * g22) * alpha * alpha * cth * sth
	- 4. * a5 * (g12 * g12) * alpha * alpha * cth2
	- 8. * a5 * g21 * g11 * g12 * g22 * alpha * alpha * cth * sth
	- 8. * a5 * g22 * g12 * alpha * alpha * cth * sth
	- 4. * a5 * (g21 * g21) * alpha * alpha * cth * sth
	- 4. * a5 * g21 * g11 * alpha * alpha * sth2
	- 4. * a5 * (g11 * g11) * (g22 * g22) * alpha * alpha * cth * sth
	- 4. * a5 * (g12 * g12) * (g21 * g21) * alpha * alpha * cth * sth
	- 4. * a5 * g11_3 * g21 * alpha * alpha * cth2
	- 4. * a5 * (g22 * g22) * alpha * alpha * cth * sth
	- 4. * a5 * (g11 * g11) * g12 * g22 * alpha * alpha * cth2
	- 8. * a5 * g21 * g11 * alpha * alpha * cth * sth
	- 4. * a5 * g22 * g12 * (g21 * g21) * alpha * alpha * sth2
	- 4. * a5 * (g12 * g12) * alpha * alpha * cth * sth
	- 4. * a5 * g21 * g11 * (g22 * g22) * alpha * alpha * sth2;
}

}

#endif // !pose_from_point_tangents_2_hxx_


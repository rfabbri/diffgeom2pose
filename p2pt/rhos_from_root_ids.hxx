#include "poly.h"
#include "pose_from_point_tangents_2_fn_t_for_root.hxx"
#include <boost/math/tools/roots.hpp>

namespace P2Pt {

template <typename T>
void
pose_poly<T>::
rhos_from_root_ids(
	const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
	T (*output)[8][ROOT_IDS_LEN] /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts, ts_len} */
)
{
	T (&rhos1)[ROOT_IDS_LEN]       = (*output)[0];
	T (&rhos1_minus)[ROOT_IDS_LEN] = (*output)[1];
	T (&rhos1_plus)[ROOT_IDS_LEN]  = (*output)[2];
	T (&rhos2)[ROOT_IDS_LEN]       = (*output)[3];
	T (&rhos2_minus)[ROOT_IDS_LEN] = (*output)[4];
	T (&rhos2_plus)[ROOT_IDS_LEN]  = (*output)[5];
	T (&ts)[ROOT_IDS_LEN]          = (*output)[6];
	T (&ts_len)                    = (*output)[7][0];

	static int ts_end;
	ts_end = 0;

	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		if (root_ids[i] == 1) {
			// TODO: implement fzero using (t_vector(i) + t_vector(i+1))/2; maybe use Newton's method later
			static T t_ref_arr[11];
			static T t_ref;

			pose_from_point_tangents_2_fn_t_for_root((t_vector[i] + t_vector[i+1])/2, &t_ref_arr);

			//using namespace boost::math::tools;
			//unsigned long max_iter = 100;
			//auto test = toms748_solve(&P2Pt::pose_poly<double>::pose_from_point_tangents_2_fn_t_for_root, t_vector[i], t_vector[i+1], eps_tolerance<T>(), max_iter);

			t_ref = t_ref_arr[0];

			//pose_from_point_tangents_2_fn_t_for_root(t_ref);

			// TODO: Possibly optimize the size of `ts[]`.
			// What is the max number of 1s than can appear in `root_ids[]`?
			ts[ts_end++] = t_ref;
		}
	}
	ts_len = ts_end;

	#if 0 // FOR TEST ONLY
	ts_end = 4;
	ts_len = 4;
	ts[0] = -0.276012891405233;
	ts[1] = -0.134802317714200;
	ts[2] = 0.481519295339129;
	ts[3] = 0.711198234568026;
	#endif

	T t_stddev = t_vector[1] - t_vector[0];

	#define ALPHA_TS_COS(x) (2 * alpha * (x) * cos(theta))
	#define ALPHA_TS_SIN(x) (-2 * alpha * (x) * sin(theta))
	#define BETA_TS_SIN(x) (beta * (1 - (x) * (x)) * sin(theta))
	#define BETA_TS_COS(x) (beta * (1 - (x) * (x)) * cos(theta))
	#define TS_DEN(x) (1 + (x) * (x))

	//% Each root is now ts(i), plus minus t_stddev.
	//% Now get rho1(t):

	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i];
		rhos1[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1[i] /= TS_DEN(ts_new);
		rhos2[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i] - t_stddev;
		rhos1_minus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_minus[i] /= TS_DEN(ts_new);
		rhos2_minus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_minus[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i] + 2 * t_stddev;
		rhos1_plus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_plus[i] /= TS_DEN(ts_new);
		rhos2_plus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_plus[i] /= TS_DEN(ts_new);
	}

}

}


